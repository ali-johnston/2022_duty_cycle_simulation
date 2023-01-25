
#########################################################
## title: functions to simulate sound distributions and duty cycles
## purpose: proof of concept
## author: ali johnston
## date: 2022.12.10


#########################################################
## script headers


# -------------------------------------------------------
# package loading
library(tidyverse)
library(forecast)



#########################################################
## create the simulation functions

# these run entirely in chunks of 'units', which are
# arbitrary time for these functions. they are only
# defined outside the functions



#########################################################
## simulation to create units of sound for different
## individuals. 

# sim_sound 
# simulates the sound distribution of a single organism across a time period
# p_sound: 			proportion of units that the organism is making identifiable sounds
# p_silence_consec: proportion of _consecutive_ units that there is silence
# number_units: 	total number of units in this sequence. must be less than 1-p_sound
# cluster_strength: log of the clustering value. 
#					the exp of this is the multiplier for relative probabilites of sound
# 					in a given unit, on the logit scale. 
#					a small value leads to very little difference on logit scale and on prob scale
#					a large value (e.g. 2) leads to bimodal probability distribution between 0 and max.  
# 					value between:
#							NA (no clustering - default)
#							0 (negligable clustering)
# 							1 (moderate clustering)
# 							2 (extreme clustering)
# 					if there is a low prevalence of sound within the available time, the clustering 
# 					is not very powerful with this mechanism. 
# random_start: 	randomise the start time when organism starts sound. if FALSE it defaults to starting with sound
# seed: 			a seed if required. 
# n_animals: 		the number of animals to simulate with the same parameters. defaults to 1. 
# plot_sound: 		if true, plots a quick diagram of the sound distributions
# order: 			the order of moving average for the peaks/troughs of any clustering in sound. 41 = 1 hour. so on average one peak and trough per hour. 
# verbose:			if TRUE, prints out information about the run

sim_sound <- function(p_sound = 0.20,
					p_silence_consec = 0.50,
					number_units = 480,
					cluster_strength = NA,
					random_start = TRUE, 
					seed = NULL, 
					n_animals = 1, 
					plot_sound = TRUE, 
					order = 41, 
					verbose = FALSE,
					...){

	# print statement about location
	if(verbose) print("---------- SIMULATING SOUNDS -----------")


	# check p_silence is not implausibly high
	if(p_silence_consec > (1-p_sound)){
		p_silence_consec <- 1 - p_sound
		print(paste("p_silence_consec adjusted to", 1-p_sound, "as it cannot be more than 1-p_sound"))
	}

	# set seed
	if(is.null(seed)) seed <- ceiling(runif(1, 1, 1000))
	set.seed(seed)

	if(verbose){
		print(paste("p_sound =", round(p_sound, 3)))
		print(paste("p_silence_consec =", round(p_silence_consec, 3)))
		print(paste("cluster_strength =", cluster_strength))
		print(paste("n_animals =", n_animals))
		print(paste("seed =", seed))
	}

	# -------------------------------------------------------
	# define (in units) number of consec silence, interspersed silence, and sound

	silence_units_consec <- floor(p_silence_consec * number_units)
	sound_units <- floor(p_sound * number_units)
	silence_units_interspersed <- number_units - silence_units_consec - sound_units
	potential_sound_units <- number_units - silence_units_consec


	# -------------------------------------------------------
	# set up mini-function to simulate a single animal
	sim_1_animal <- function(dummy = 1, potential_sound_units, sound_units, cluster_strength = 0, order = 41, seed = 1, ...){

		set.seed(dummy + (seed*1000))

		if(sound_units == potential_sound_units){
			sound_locs_1n <- rep(1, sound_units)
		}

		if(sound_units < potential_sound_units){

			# clustered probability. 
			# parameter cluster_strength defines the degree of clustering
			# 'gaussian random field' in temporal dimension
			# calculate moving average of random values
			# then standardise so max and min are 0 and 1. 

			rn <- rnorm(potential_sound_units + order - 1, 0, 1) |>
					forecast::ma(order = order)
			rn <- rn[!is.na(rn)]
			if((max(rn) - min(rn)) == 0) rn <- rep(1, potential_sound_units)
			if((max(rn) - min(rn)) > 0) rn <- (rn-min(rn))/(max(rn)-min(rn))

			# adjust the probabilities by multiplying the standardised moving average by the 
			# clustering strength. 

			# if cluster_strength = 0, they are all equal probability (random distribution of sound)
			# if cluster_strength = 1, sound is clustered, with probabilities going to 0 for some units of time

			# Q: do we want a distribution that is more even than random? 
			# Q: what patterns do real species show?
			if(is.na(cluster_strength)) logit_probs <- 0 + (rn-0.5)
			if(!is.na(cluster_strength)) logit_probs <- 0 + (rn-0.5) * exp(cluster_strength)
			unscaled_probs <- exp(logit_probs) / (1 + exp(logit_probs))
			m_probs <- unscaled_probs / sum(unscaled_probs)

			sound_loc_indices <- sample(1:potential_sound_units, size = sound_units, replace = FALSE, prob = m_probs) |> sort()
			sound_locs_1n <- rep(0, potential_sound_units)
			sound_locs_1n[sound_loc_indices] <- 1
		}
		return(sound_locs_1n)
	}

	# -------------------------------------------------------
	# do the simulations for all animals and combine together in 
	# a matrix with the consecutive silence
	sound_locs_list <- map(1:n_animals, sim_1_animal, 
		potential_sound_units = potential_sound_units, 
		sound_units = sound_units, 
		cluster_strength = cluster_strength, 
		seed = seed, 
		order = order)
	sound_locs <- do.call(cbind, sound_locs_list)

	sound_24 <- rbind(sound_locs, matrix(0, ncol=n_animals, nrow = silence_units_consec))


	# -------------------------------------------------------
	# if required, randomise the start time for each individual
	if(random_start){
		start_index <- sample(1:number_units, size = n_animals, replace = TRUE)
		start_n <- number_units - start_index + 1
		sound_24_newstart <- sound_24
		# loop through animals and reorder the start to be at the random start unit
		for(i in 1:n_animals){
			if(start_index[i] > 1) {
				sound_24_newstart[1:start_n[i],i] <- sound_24[start_index[i]:number_units,i]
				sound_24_newstart[(start_n[i]+1):number_units,i] <- sound_24[1:(start_index[i]-1),i]
			}
		}
		sound_24 <- sound_24_newstart
	}

	if(plot_sound) plot_sound_fn(sound_24, ...)

	return(list(sounds = sound_24, seed = seed))
}


#########################################################
## plot the simulated sounds. 

# takes the output from sim_sound

plot_sound_fn <- function(sound_24, col_sound = "darkblue", ...){
	image(sound_24, 
			yaxt = "n", ylab = "individuals", 
			xaxt = "n", xlab = "sound units (colour = sound)",
			col = c("white", col_sound), ...)
}



#########################################################
## simulate the duty cycles with input


# sim_duty
# simulates a duty cycle across a time period

# prop_recording:	the proportion of units where recording is happening
# units_per_chunk:	the number of time chunks in each recording (on) block
# number_units: 	total number of units in this sequence. must be less than 1-p_sound
# random_start: 	randomise the start time in the duty cycle. if FALSE it defaults to starting with "on"
# seed: 			a seed if required. 
# plot_duty: 		if TRUE, plots a quick diagram of the duty cycle
# plot_duty_add:	if TRUE, adds the duty cycle onto any existing plot (designed to be on a sound plot)
# verbose:			if TRUE, prints out information about the run

sim_duty <- function(prop_recording = 0.5,
					units_per_chunk = 40,
					number_units = 480,
					random_start = TRUE, 
					seed = NULL,
					plot_duty = TRUE,
					plot_duty_add = TRUE, 
					verbose = FALSE,
					...){

	# print statement about location
	if(verbose) print("---------- SIMULATING DUTY CYCLE -----------")

	if(is.null(seed)) seed <- ceiling(runif(1, 1, 1000))
	set.seed(seed)

	prop_on <- prop_recording
	units_on <- units_per_chunk
	prop_off <- 1 - prop_recording
	units_off <- round(units_per_chunk / prop_recording * prop_off, 0)


	if(verbose){
		print(paste("prop_recording =", round(prop_recording, 3)))
		print(paste("units_per_chunk =", units_per_chunk))
		print(paste("seed =", seed))
	}

	# start with defining a single on/off cycle
	single_cycle <- c(rep(1, units_on), rep(0, units_off))

	# repeat up to the total number of units (usually a 24-hour period)
	rep_no <- ceiling(number_units / length(single_cycle))
	multi_cycle <- single_cycle
	for(i in 1:rep_no) multi_cycle <- c(multi_cycle, single_cycle)

	# double up to 2 x number_units, so that you can start anywhere for random
	multi_cycle <- c(multi_cycle, multi_cycle)

	# if not randomising, use the first units
	if(!random_start) loop_day <- multi_cycle[1:number_units] 

	# if random start, pick a random starting point in the cycle
	if(random_start) {
		start <- ceiling(runif(1, 0.5, number_units + 0.5))
		loop_day <- multi_cycle[start:(start + number_units -1)]

	}

	if(plot_duty) plot_duty_fn(loop_day, n_animals = 2, add = plot_duty_add, ...)

	return(list(duty = loop_day, seed = seed))

}


#########################################################
## plot the duty cycles. 

# this is constructed to fit with the sound plots. 

plot_duty_fn <- function(loop_day, n_animals, add = TRUE, col_on = alpha("firebrick", 0.3), ...){
	duty_n_animals <- loop_day
	if(n_animals>1){
		for(i in 2:n_animals) duty_n_animals <- rbind(duty_n_animals, loop_day)
	}
	image(t(duty_n_animals), 
			yaxt = "n", ylab = "", 
			xaxt = "n", xlab = "",
			col = c(alpha("white", 0), col_on), 
			add = add)
}


#########################################################
## overall function to combine simulations of sound, duty cycles, and plots

# sim_sound_duty
# simulates the sound distribution of a single organism across a time period
# simulates a duty cycle
# combines these together and reports which animals were detected

# p_sound: 			proportion of units that the organism is making identifiable sounds
# p_silence_consec: proportion of _consecutive_ units that there is silence
# number_units: 	total number of units in this sequence. must be less than 1-p_sound
# cluster_strength: value between 0 (no clustering) to 1 (extreme clustering)
# random_start_sound:randomise the start time when organism starts sound. if FALSE it defaults to starting with sound
# seed: 			a seed if required. 
# n_animals: 		the number of animals to simulate with the same parameters. defaults to 1. 
# order: 			the order of moving average for the peaks/troughs of any clustering in sound. 41 = 1 hour. so on average one peak and trough per hour. 

# prop_recording:	the proportion of units where recording is happening
# units_per_chunk:	the number of time chunks in each recording (on) block
# random_start_duty:randomise the start time in the duty cycle. if FALSE it defaults to starting with "on"

# plot_sound_duty: 	if TRUE, plots a quick diagram of the sound distributions within the sound function
# verbose:			if TRUE, prints out information about the run

sim_sound_duty <- function(p_sound = 0.20,
					p_silence_consec = 0.50,
					number_units = 480,
					cluster_strength = 0,
					random_start_sound = TRUE, 
					seed = NULL, 
					n_animals = 1, 
					order = 41,
					prop_recording = 0.5,
					units_per_chunk = 40,
					random_start_duty = TRUE,
					plot_sound_duty = FALSE, 
					verbose = FALSE, 
					...){

sounds <- sim_sound(p_sound = p_sound,
					p_silence_consec = p_silence_consec,
					number_units = number_units,
					cluster_strength = cluster_strength,
					random_start = random_start_sound, 
					seed = seed, 
					n_animals = n_animals, 
					plot_sound = FALSE, 
					order = order, 
					verbose = verbose,
					...)

duty <- sim_duty(prop_recording = prop_recording,
					units_per_chunk = units_per_chunk,
					number_units = number_units,
					random_start = random_start_duty, 
					seed = seed, 
					plot_duty = FALSE, 
					verbose = verbose,
					...)


# -------------------------------------------------------
# plot_both
if(plot_sound_duty){
	plot_sound_fn(sounds$sounds, ...)
	plot_duty_fn(duty$duty, n_animals = n_animals, ...)
}

# -------------------------------------------------------
# combine these

no_detections <- as.matrix(t(sounds$sounds)) %*% as.matrix(duty$duty)

params = list(p_sound = p_sound, p_silence_consec = p_silence_consec, number_units = number_units,
			cluster_strength = cluster_strength, random_start_sound = random_start_sound, 
			seed_sound = sounds$seed, seed_duty = duty$seed,
			n_animals = n_animals, order = order,
			prop_recording = prop_recording, units_per_chunk = units_per_chunk,
			random_start_duty = random_start_duty)
return_list <- list(no_detections = no_detections, sounds = sounds$sounds, duty = duty$duty,
		params = params)

return(return_list)
}










#########################################################
## title: run several simulations
## purpose: proof of concept for duty cycle plan
## author: ali johnston
## date: 2022.12.11


#########################################################
## script headers


# -------------------------------------------------------
# package loading
library(tidyverse)
library(RColorBrewer)

# -------------------------------------------------------
# folder paths

wd <- "/Users/Alison/Documents/REPOS/2022_duty_cycle_simulation/"
setwd(wd)


#########################################################
## load in functions

source(paste0(wd, "code/step1_simulation_functions.R"))


#########################################################
## test simulation functions for single runs

# add clustering of sounds
test_sound2 <- sim_sound(p_sound = 0.20,		
					p_silence_consec = 0.20, 	
					number_units = 480, 		
					cluster_strength = 1,		# create clustering of sounds within sound window
					random_start = FALSE, 		
					n_animals = 5, 
					order = 61)

test_sound1 <- sim_sound(p_sound = 0.20,		# 20% of chunks with sound occuring
					p_silence_consec = 0.50, 	# consecutive 50% of each day with no sound
					number_units = 480, 		# 3 minute chunks. 480 over 24-hour period
					cluster_strength = 0,
					random_start = FALSE, 		# don't randomise start (easier to see what it's doing)
					n_animals = 5, 
					order = 61)

# add clustering of sounds
test_sound2 <- sim_sound(p_sound = 0.20,		
					p_silence_consec = 0.80, 	
					number_units = 480, 		
					cluster_strength = 1,		# create clustering of sounds within sound window
					random_start = FALSE, 		
					n_animals = 5, 
					order = 61)

# randomise start time of sound
test_sound3 <- sim_sound(p_sound = 0.20,		
					p_silence_consec = 0.50, 	
					number_units = 480, 		
					cluster_strength = 1,		
					random_start = TRUE, 		# randomise when individuals make sound
					n_animals = 10)

# simulate a duty cycle
test_duty <- sim_duty(prop_recording = 0.20,	# recording for 20% of the day
					units_per_chunk = 40,		# recording in chunks of 40 units (2 hours)
					number_units = 480,
					random_start = TRUE) 		# randomise start time


par(mfrow = c(2,1))
# simulate sound and duty cycles together
test_sound_duty1 <- sim_sound_duty(p_sound = 0.10,
					p_silence_consec = 0.70,
					number_units = 480,
					cluster_strength = 1,
					order = 41,
					random_start_sound = TRUE, 
					n_animals = 1000, 
					prop_recording = 0.20,
					units_per_chunk = 40,
					random_start_duty = TRUE,
					plot_sound_duty = FALSE)

test_sound_duty2 <- sim_sound_duty(p_sound = 0.10,
					p_silence_consec = 0.70,
					number_units = 480,
					cluster_strength = 0,
					random_start_sound = TRUE, 
					n_animals = 1000, 
					prop_recording = 0.20,
					units_per_chunk = 40,
					random_start_duty = TRUE,
					plot_sound_duty = FALSE)


# use seed to get the same simulations:
par(mfrow = c(2,1))
# simulate sound and duty cycles
test_sound_duty1 <- sim_sound_duty(p_sound = 0.20,
					p_silence_consec = 0.50,
					number_units = 480,
					cluster_strength = 0,
					random_start_sound = TRUE, 
					n_animals = 20, 
					prop_recording = 0.20,
					units_per_chunk = 40,
					random_start_duty = TRUE,
					plot_sound_duty = TRUE, 
					seed = 1234)

test_sound_duty2 <- sim_sound_duty(p_sound = 0.20,
					p_silence_consec = 0.50,
					number_units = 480,
					cluster_strength = 0,
					random_start_sound = TRUE, 
					n_animals = 20, 
					prop_recording = 0.20,
					units_per_chunk = 40,
					random_start_duty = TRUE,
					plot_sound_duty = TRUE,
					seed = 1234)

# use verbose to get a report of parameters used for each run
test_sound_duty3 <- sim_sound_duty(verbose = TRUE)


# output of sim_sound_duty is a list with: 
# no_detections: 	which individuals have been detected and how many times
# sounds: 			the output from the sound simulations
# duty: 			the output from the duty simulations
# params:			a list of parameters used

names(test_sound_duty2)

# to get the proportion of individuals detected: 
length(test_sound_duty2$no_detections) / test_sound_duty2$params$n_animals










#########################################################
## set up simulation parameters for varying characteristics
## of sound distribution across a 24-hour period


# for now the simulations are in units of 3 mins

# set up variation in proportion of the day with sound
single_unit_mins <- 3
tot_mins <- 24*60
tot_units <- tot_mins / single_unit_mins
units_sound <- seq(1, round(tot_units*0.25), by = 5)

# set up variation in how many units of consecutive silence
units_silence_consec <- seq(0.5, 0.95, by = 0.05)*tot_units |> round()
units_sound_contained <- tot_units - units_silence_consec

# set up distribution in gaps between units with sound



#########################################################
## set duty cycle parameters 

## for now just picking one as a test case
## eventually this will be automated in an array of 
## simulations

propn_sampling <- 0.25
number_units_per_chunk <- 40 
# with 3-minute units:
# 20 units = 1 hour



#########################################################
## set up results matrix and loop through simulations

res_mat <- matrix(NA, nrow = length(units_sound), ncol = length(units_silence_consec))

for(kkk in 1:length(units_sound)){

	viable_gaps <- (tot_units - units_silence_consec) >= units_sound[kkk]
	if(any(viable_gaps)){
		for(jjj in which(viable_gaps)){

			n_animals <- 100

			run_sim <- sim_sound_duty(p_sound = units_sound[kkk] / tot_units,
								p_silence_consec = units_silence_consec[jjj] / tot_units,
								number_units = tot_units,
								cluster_strength = 1,
								seed = (kkk*1000) + jjj, 
								n_animals = n_animals, 
								order = 41,
								prop_recording = propn_sampling,
								units_per_chunk = number_units_per_chunk)

			detected <- run_sim$no_detections > 0
			res_mat[kkk, jjj] <- as.numeric(sum(detected) / length(detected))

		}
	}
}


res_long <- data.frame(units_sound = rep(units_sound, length(units_silence_consec)),
						units_silence_consec = sort(rep(units_silence_consec, length(units_sound))),
						prop_detected = as.vector(res_mat)) |>
				mutate(prop_sound = units_sound / tot_units) |>
				mutate(prop_sound_contained = 1 - (units_silence_consec / tot_units))


ggplot(data=res_long, aes(x = prop_sound, y = prop_sound_contained)) + 	
	geom_tile(aes(fill = prop_detected)) + 
	scale_fill_gradientn(limits = c(0, 1), colours = brewer.pal(8, "YlGnBu")) + 
	xlab("Proportion of day with sound") + 
	ylab("Proportion of day containing all sound") + 
	theme_bw() + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())





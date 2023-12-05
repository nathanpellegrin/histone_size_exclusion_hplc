setwd("~/Documents/BTECH 609 Histone Size Exclusion Lab")

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(r2symbols)

# Assuming your data is in a CSV file named "data.csv"
data <- read.csv("HPLC_histones_size_exclusion.csv")

# Reshape the data to long format
data_long <- data %>%
  gather(key = "variable", value = "value", -Time) %>%
  mutate(
    Sample = substr(variable, 1, nchar(variable) - 5),
    Wavelength = paste0(substr(variable, nchar(variable) - 3, nchar(variable)), "Å"),
    Wavelength = sub("A", "", Wavelength),
    value = as.numeric(value)
  ) %>%
  select(-variable)

data_long[is.na(data_long$value),]

sample_dat <- sample_n(data_long, 100)

write.csv(sample_dat, "sample.csv", row.names = FALSE)

Greek_Letters <- data.frame(
Name= c("Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Zeta", "Eta", "Theta", "Iota", "Kappa", "Lambda", "Mu", "Nu", "Xi", "Omicron", "Pi", "Rho", "Sigma", "Tau", "Upsilon", "Phi", "Chi", "Psi", "Omega"),
Lowercase=c("α", "β", "γ", "δ", "ε", "ζ", "η", "θ", "ι", "κ", "λ", "μ", "ν", "ξ", "ο", "π", "ρ", "σ", "τ", "υ", "φ", "χ", "ψ", "ω"),
Uppercase=c("Α", "Β", "Γ", "Δ", "Ε", "Ζ", "Η", "Θ", "Ι", "Κ", "Λ", "Μ", "Ν", "Ξ", "Ο", "Π", "Ρ", "Σ", "Τ", "Υ", "Φ", "Χ", "Ψ", "Ω"))

get_gl <- function(name, case = "lower") {
  Greek_Letters <- data.frame(
    Name = c("Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Zeta", "Eta", "Theta", "Iota", "Kappa", "Lambda", "Mu", "Nu", "Xi", "Omicron", "Pi", "Rho", "Sigma", "Tau", "Upsilon", "Phi", "Chi", "Psi", "Omega"),
    Lowercase = c("α", "β", "γ", "δ", "ε", "ζ", "η", "θ", "ι", "κ", "λ", "μ", "ν", "ξ", "ο", "π", "ρ", "σ", "τ", "υ", "φ", "χ", "ψ", "ω"),
    Uppercase = c("Α", "Β", "Γ", "Δ", "Ε", "Ζ", "Η", "Θ", "Ι", "Κ", "Λ", "Μ", "Ν", "Ξ", "Ο", "Π", "Ρ", "Σ", "Τ", "Υ", "Φ", "Χ", "Ψ", "Ω")
  )
  
  # Find the closest match
  matches <- grepl(name, Greek_Letters$Name, ignore.case = TRUE)
  
  if (any(matches)) {
    match_index <- which(matches)[1] # Taking the first match
    if (tolower(case) == "lower") {
      return(Greek_Letters$Lowercase[match_index])
    } else {
      return(Greek_Letters$Uppercase[match_index])
    }
  } else {
    return(NA) # Return NA if no match is found
  }
}


ggplot(data_long, aes(x = Time, y = value, color = Wavelength)) +
  geom_line() + # Add line plot
  #facet_wrap(~ Sample) + # Facet by Sample
  facet_wrap(~factor(Sample, levels=c('SEC_Standard','Histone','Histone_Heat_Denatured',
                                      'GuHCl_Denature','GuHCL_Refolded')),ncol=2)
  theme_minimal() + # Optional: Use a minimal theme
  labs(title = "Value over Time by Sample and Wavelength",
       x = "Time",
       y = "Value",
       color = "Wavelength") +
  theme(legend.position = "bottom") # Adjust legend position

#  https://stackoverflow.com/questions/66053185/dual-y-axis-while-using-facet-wrap-in-ggplot-with-varying-y-axis-scale



data_wv <- data_long %>%
  mutate(Wavelength = gsub("Å", "A", Wavelength)) %>%
  spread(key = Wavelength, value = value) 
colnames(data_wv) <- c("Time","Sample","w214","w280")

# Function factory for secondary axis transforms
train_sec <- function(from, to) {
  from <- range(from)
  to   <- range(to)
  # Forward transform for the data
  forward <- function(x) {
    rescale(x, from = from, to = to)
  }
  # Reverse transform for the secondary axis
  reverse <- function(x) {
    rescale(x, from = to, to = from)
  }
  list(fwd = forward, rev = reverse)
}

# Learn the `from` and `to` parameters
sec <- train_sec(data_wv$w280, data_wv$w214)


ggplot(data_wv, aes(x = Time)) +
  geom_line(aes(y=w214), col="blue") + # Add line plot
  geom_line(aes(y=sec$fwd(w280)), col="red") +
  scale_y_continuous(name="wavelength", sec.axis=sec_axis(~sec$rev(.), name="A280")) +
  facet_wrap(~ Sample) + # Facet by Sample
  #facet_wrap(~factor(Sample, levels=c('SEC_Standard','Histone','Histone_Heat_Denatured',
  #                                    'GuHCl_Denature','GuHCL_Refolded')))  
  theme_minimal() + # Optional: Use a minimal theme
  labs(title = "Value over Time by Sample and Wavelength",
       x = "Time",
       y = "Wavelength") +
  theme(legend.position = "bottom") # Adjust legend position


# get max intensities for each sample
max_wl <- data_long %>%
  group_by(Sample, Wavelength) %>%
  summarize(MaxValue = max(value), MinValue=min(value)) %>%
  pivot_wider(id_cols = Sample, names_from=Wavelength, values_from = c(MaxValue,MinValue) ) 
colnames(max_wl) <- c("Sample","w214_max","w280_max","w214_min","w280_min")
max_wl$scale_factor <- max_wl$w214_max/max_wl$w280_max

# Join the max values with the original dataset and scale 214Å values
data_wv <- data_wv %>%
  left_join(max_wl, by = "Sample") %>%
  mutate(scale2=ifelse(Sample=="GuHCl_Denature",.25,.97)) %>%
  mutate(w280_scaled = w280*scale_factor*scale2) %>% 
  select(-scale_factor, -w214_max, -w280_max, -w214_min, -w280_min)

data_wv2 <- data_wv %>% filter(Sample=='Histone', Time<7.5 & Time >2.5) %>%
        mutate(Sample='Legend', w214=0, w280_scaled=50) %>% 
        bind_rows(data_wv)

max_wl <- max_wl %>% mutate(lbl=paste("w280 scaler:",round(scale_factor,3)))  %>%
        bind_rows(data.frame(Sample='Legend',lbl=" ",lbl14="214Å", lbl80="280Å"))

ggplot(data_wv2, aes(x = Time)) +
  geom_line(aes(y=w214), col="blue") + # Add line plot
  geom_line(aes(y=w280_scaled), col="red") +
  scale_y_continuous(limits=c(-50,105)) +
  facet_wrap(~factor(Sample, levels=c('SEC_Standard','Legend','Histone','Histone_Heat_Denatured',
                                      'GuHCl_Denature','GuHCl_Refolded')), drop=FALSE, ncol=2) +                                                       # Add individual text to plot
  geom_text(data = max_wl, size=3, mapping = aes(x = 9, y = 50, label = lbl80)) +
  geom_text(data = max_wl, size=3, mapping = aes(x = 9, y = 0, label = lbl14)) +
  geom_text(data = max_wl, size=3, mapping = aes(x = 2.5, y = 90, label = lbl)) +
  theme_minimal() + # Optional: Use a minimal theme
  labs(title = "Intensity of Absorption over Time by Sample and Wavelength",
       x = "Time (minutes)",
       y = "Absorption Intensity") 



# this function regenerates ggplot2 default colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

color_lty_cross = expand.grid(
  ltypes = 1:6,
  colors = gg_color_hue(10)
  ,stringsAsFactors = F)

ggplot(data_long, aes(x = Time, y = value, color = Sample, lty=Sample)) +
  geom_line(size=.75) + 
  facet_wrap(~Wavelength, ncol=1, scales="free") +
  #scale_color_manual(values = color_lty_cross$colors[1:5]) +
  scale_linetype_manual(values = color_lty_cross$ltypes[1:5]) +
theme_minimal() + 
  labs(title = "Absorption Intensity over Time by Wavelength and Sample",
       x = "Time (minutes)",
       y = "Intesity",
       color = "Sample: ") +
  theme(legend.position = "right") # Adjust legend position


sec_std <- data.frame(
molecule=c('Thyroglobulin',paste0(get_gl("gamma"),'-Globulin'),'Albumin','Riboncls-A',paste0(get_gl("rho"),'-ABA')),
weight=c(670000,150000,44300,13700,137),
wgt_kd=c('670 kDa', '150 kDa', '44.3 kDa','13.7 kDa','0.137 kDa'),
elution_time=c(6.8,10.9,12.8,14.2,17.0) )




ggplot(sec_std, aes(x=elution_time, y=log(weight))) +
    geom_point(aes(color=molecule), size=5) + geom_line() + #geom_smooth(method="loess") +
    geom_text(mapping = aes(label = wgt_kd), nudge_x=1) +
    theme_minimal() +
    labs(title = "Log Molecular Weight (Da) and Elution Time of SEC Standards",
     x = "Elution time (minutes)",
     y = "log(g/l)",
     color = "Molecule: ")

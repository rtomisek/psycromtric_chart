#! /usr/bin/Rscript

#
# Copyright (C) 2025 by Randall Tomisek  <rtomisek@gmail.com>
#
# Permission to use, copy, modify, and/or distribute this software 
# for any purpose with or without fee is hereby granted.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL
# WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL
# THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
# OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
# NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION
# WITH THE USE OR PERFORMANCE OF THIS SOFTWARE
#

# An R script to generate psychrometric charts


# A function to return saturated water vapor pressure at temperature t
p_sat <- function( t) {

   if( t < 213.15 )
   {
      return( -1 )
    } else if( t < 273.15 )
   {
      A <- -0.7297593707E-5
      B <-  0.5397420727E-2
      C <-  0.2069880620E2
      D <- -0.6042275128E4
   } else if( t < 322.15 )
   {
      A <-  0.1255001965E-4
      B <- -0.1923595289E-1
      C <-  0.2705101899E2
      D <- -0.6344011577E4
   } else if( t < 373.15 )
   {
      A <-  0.1246732157E-4
      B <- -0.1915465806E-1
      C <-  0.2702388315E2
      D <- -0.6340941639E4
   } else if( t < 423.15 )
   {
      A <-  0.1204507646E-4
      B <- -0.1866650553E-1
      C <-  0.2683629403E2
      D <- -0.6316972063E4
   }  else if( t <= 473.15 )
   {
      A <-  0.1069730183E-4
      B <- -0.1698965754E-1
      C <-  0.2614073298E2
      D <- -0.6220781230E4
   } else
   {
      return( -1)
   }
    
   a <- ( (A*t) + B )*t + C + D/t
   
   return( 1000 * exp( a ) )
}

# used to gererate data for rh lines
Wrh <- function( t, rh, P ) {
    return( 0.62198*p_sat(t+273.15)*rh / (P - rh*p_sat(t+273.15)) )
}
vWrh <- Vectorize(Wrh, "t")

# used to generate enthalpy lines
We <- function( t, H, P) {
    return( (H - 1.006*t)/(2501 + 1.805*t) )
}
vWe <- Vectorize(We, "t")


# parameters this chart
P <- 101325  # pascals, pressure for chart, usually 1 atm

Tmin <- -5.0
Tmax <- 55.0

Wmin <- 0.0
Wmax <- 0.032

title1 <- "Psychrometric Chart"                                       # make strings for main label
title2 <- sprintf("Pressure %s pascals", format(P, trim=TRUE, big.mark=",") )


library(ggplot2)

p <- ggplot( )  +
     annotate( "text", x=5, y=0.030, label=title1, size=5) +                                       # main title
     annotate( "text", x=5, y=0.029, label=title2, size=3) +
     theme( panel.grid = element_line( color='black', linetype=2, linewidth=0.05 ),
            panel.background = element_rect(fill = 'white', colour='black'  )) +
     xlab("Dry Bulb Temperature, °C") +
     scale_x_continuous(n.breaks=21, minor_breaks = scales::breaks_width(1), limits=c(Tmin, Tmax) ) +
     ylab("Humidity Ratio, W") +
     scale_y_continuous(n.breaks=18, limits=c(Wmin, Wmax), position="right") +
     geom_function( fun = function(x) vWrh( x, 1.0, P) ) +                                  # rh = 1.0 it the saturation line
     geom_function( fun = function(x) vWrh( x, 0.8, P), linewidth=0.1 ) +                   # is should be plotted bolder that the others
     geom_function( fun = function(x) vWrh( x, 0.6, P), linewidth=0.1 ) +
     geom_function( fun = function(x) vWrh( x, 0.4, P), linewidth=0.1 ) +                   # add the whateve rh lines you want
     geom_function( fun = function(x) vWrh( x, 0.2, P), linewidth=0.1 ) +
     geom_function( fun = function(x) vWrh( x, 0.1, P), linewidth=0.1 ) +
     annotate( "text", x=28, y=0.0185, label="rh=80%", size=2, fontface='italic', angle=60) +
     annotate( "text", x=31.8, y=0.0172, label="rh=60%", size=2, fontface='italic', angle=60) +  # labels for rh lines
     annotate( "text", x=36.5, y=0.015, label="rh=40%", size=2, fontface='italic', angle=55) +
     annotate( "text", x=45.2, y=0.0117, label="rh=20%", size=2, fontface='italic', angle=50) +
     annotate( "text", x=52.55, y=0.0084, label="rh=10%", size=2, fontface='italic', angle=36) +
     geom_function( fun = function(x) vWe( x, 10, P), xlim=c(-4.0,55.0), linewidth=0.1 ) +
     geom_function( fun = function(x) vWe( x, 20, P), xlim=c(1.2,55.0), linewidth=0.1 ) +
     geom_function( fun = function(x) vWe( x, 30, P), xlim=c(7.0,55.0), linewidth=0.1 ) +
     geom_function( fun = function(x) vWe( x, 40, P), xlim=c(12.0,55.0), linewidth=0.1 ) +
     geom_function( fun = function(x) vWe( x, 50, P), xlim=c(15.0,55.0), linewidth=0.1 ) +      # the entalpy lines
     geom_function( fun = function(x) vWe( x, 60, P), xlim=c(17.5,55.0), linewidth=0.1 ) +
     geom_function( fun = function(x) vWe( x, 70, P), xlim=c(20.0,55.0), linewidth=0.1 ) +
     geom_function( fun = function(x) vWe( x, 80, P), xlim=c(22.5,55.0), linewidth=0.1 ) +
     geom_function( fun = function(x) vWe( x, 90, P), xlim=c(25.0,55.0), linewidth=0.1 ) +
     geom_function( fun = function(x) vWe( x, 100, P), xlim=c(27.5,55.0), linewidth=0.1 ) +
     geom_function( fun = function(x) vWe( x, 110, P), xlim=c(29.5,55.0), linewidth=0.1 ) +
     geom_function( fun = function(x) vWe( x, 120, P), xlim=c(31.0,55.0), linewidth=0.1 ) +
     geom_function( fun = function(x) vWe( x, 130, P), xlim=c(32.0,55.0), linewidth=0.1 ) +
     annotate( "text", x=15, y=0.02, label="Enthalpy kJ/kg", size=3.3, angle=60 ) +
     annotate( "text", x=-4.0, y=0.0055, label="10", size=2.5) +
     annotate( "text", x=0.95, y=0.0077, label="20", size=2.5) +
     annotate( "text", x=6.5, y=0.0095, label="30", size=2.5) +
     annotate( "text", x=11.6, y=0.0113, label="40", size=2.5) +
     annotate( "text", x=14.8, y=0.014, label="50", size=2.5) +
     annotate( "text", x=17.3, y=0.0169, label="60", size=2.5) +
     annotate( "text", x=20.0, y=0.0197, label="70", size=2.5) +                                 # and the enthapy labels
     annotate( "text", x=22.5, y=0.0226, label="80", size=2.5) +
     annotate( "text", x=24.8, y=0.0255, label="90", size=2.5) +
     annotate( "text", x=27.3, y=0.0285, label="100", size=2.5) +
     annotate( "text", x=29.5, y=0.0315, label="110", size=2.5) +
     annotate( "text", x=38.0, y=0.032, label="120", size=2.5) +
     annotate( "text", x=47.0, y=0.032, label="130", size=2.5) +
     # An equation for the user to calculate specific volume
     annotate( "text", x=-3, y=0.0175, label="Specific volume given by", size=2.5) +
     annotate( "text", x=-3, y=0.0165, label=expression(V[sp]==frac(0.287055*T,P)~(1+1.6078*W) ), size=2.5)


ggsave("psy.pdf", plot=p, units="in", width=10, height=8)

  

#! /usr/bin/Rscript
#
# Copyright (C) 2025 by [copyright holder] <[email]>
#
# Permission to use, copy, modify, and/or distribute this software for any purpose with or without fee is hereby granted.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
# SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE,
# DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
# ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE
#

# an R script to generate psychrometric charts

################################################################################
#    EDIT TO YOUR PARAMETERS FOR THE CHART                                     #

P <- 101325                    # pascals, pressure for chart, usually 1 atm

Tmin <- -5.0                   # degrees C
Tmax <- 55.0

Wmin <- 0.0                    
cut.frac <- 3

rh.values <- c(0.8, 0.6, 0.4, 0.2, 0.1)    # Enter rh fractions you want plotted
                                           # The fraction 1.0 is aways plotted automatically
VspNote <- TRUE


Tstep.l <- 5.0                 # 'x' axis for labels 
Tstep.g <- 1.0                 # 'x' for grid
nWstep.l <- 10                 # for 'y' axis labels
nWstep.g <- 20
nHsteps <- 20

lw <- 0.01                                 # for most curves

#    End of user parameters                                                    #
################################################################################

################################################################################
#               Psychrometric functions                                        #

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

# used to gererate rh lines
Wrh <- function( t, rh, P ) {
    return( 0.62198*p_sat(t+273.15)*rh / (P - rh*p_sat(t+273.15)) )
}
vWrh <- Vectorize(Wrh, vectorize.args="t")           # must use vectorize in loops

# used to generate enthalpy lines
We <- function( t, H, P) {
    return( (H - 1.006*t)/(2501 + 1.805*t) )
}
vWe <- Vectorize( We, vectorize.args="t")

h <- function( t, W ) {
   return( 1.005*t + W*( 2501 + 1.805*t) ) 
}
#                   End psychrometric functions                                #
################################################################################

# preliminary calculations and useful function definitions

MaxW <- Wrh( Tmax, 1.0, P )                    # Estimate Wmax
Wmax <- MaxW/cut.frac                          # 
Wmax <- as.integer(Wmax*10000+0.5)/10000       # round off to 4 decimal places

Wstep.l <- as.integer(10000*Wmax/nWstep.l)/10000
Wstep.g <- Wmax/nWstep.g

sr <- (Tmax-Tmin)/(Wmax-Wmin)                  # coordinate correction for slopes

rh.x.offset <- (Tmax-Tmin)/500                 # Offset labels from lines. AVs
rh.y.offset <- (Wmax-Wmin)/200

rh.ref <- function( t ) {                                # A (not plotted) reference line to position rh labels
    return( Wmax + ((Wmin-Wmax)/(Tmax-Tmin))*(t-Tmin) )
}

Hmin <- h(Tmin, Wmin);     Hmax <- h(Tmax, Wmax)                          # Calculate the enthapy lines to
Hsteps <- round( (((Hmax-Hmin)/nHsteps)+2)/5 )*5                          #
Hmin <- Hsteps*round(Hmin/Hsteps);     Hmax <- Hsteps*round(Hmax/Hsteps)  # be plotted
H.values <- trunc(seq( Hmin+Hsteps, Hmax, Hsteps ) )                      #

pt <- uniroot( function(x) Wmax-Wrh(x, 1.0, P), c(Tmin, Tmax), tol=0.001)
r1 <- c( Tmin, Wrh(Tmin, 1.0, P)+0.001 )
r2 <- c( pt[[1]]-1, Wrh(pt[[1]], 1.0, P) )
slope <- (r2[2]-r1[2])/(r2[1]-r1[1])
H.ref <- function( t ) {                                          # A (not plotted) reference line to position H labels
      return( r1[2] + (slope)*(t-r1[1]) ) 
}


# Start chart commands
title.s <- sprintf("Psychrometric Chart\nPressure %s pascals", format(P, big.mark=",") )    # Make a string for the  main label

pdf( "psy.pdf", width=10, height=8.0 )

par( mar=c(5,2,2,5) )
curve( vWrh( x , 1.0, P), lwd=3,                         # First plot the saturation line, rh=1.0
       xlim=c(Tmin, Tmax), ylim=c(Wmin, Wmax),
       xaxt='n', yaxt='n',
       xlab="", ylab="" )

axis(side=1, at=seq(Tmin, Tmax, Tstep.l), cex.axis=0.8 )              # 'x' axis
mtext("Temperature \u00B0C", side=1, line=3 )

axis(side=4, at=seq(Wmin, Wmax, Wstep.l ), las=2, cex.axis=0.75 )     # 'y' axis and label on right
mtext("Humidity Ratio, W", side=4, line=3)

abline( v=seq(Tmin, Tmax, Tstep.g), lwd=lw, col="light gray" )        # The grid
abline( h=seq(Wmin, Wmax, Wstep.g), lwd=lw, col="light gray" )        #

text( x=Tmin+Tstep.g, y=Wmax-Wstep.g, label=title.s, pos=4 )          # The main label

for( t in seq(Tmin+Tstep.l, Tmax-Tstep.l, Tstep.l) ) {
     text( x=t, y=vWrh(t, 1.0, P)+Wstep.g/2, label=format(t), cex=0.5 )
}                                                                     # 'x' axis Temperature cues on saturation line

for( r in rh.values ) {
    curve( vWrh(x, r, P), lwd=lw, add=TRUE )
    pt <- uniroot( function(x) rh.ref(x)-vWrh(x, r, P), c(Tmin,Tmax), tol=0.001 )
    s <- ( vWrh(pt[[1]]-5.0, r, P)-vWrh(pt[[1]]+0.0, r, P))/5         # Local slope of line. AV 5
    a <- -atan(s*sr)*180/pi                                           # The angle for the label
    lb <- paste0( "rh=", format(r*100), "%" )
    text( x=pt[[1]]+rh.x.offset, y=vWrh(pt[[1]], r, P)-rh.y.offset, label=lb, cex=0.5, font=3, srt=a )
}

for( h in H.values ) {
    pt <- uniroot( function(x) H.ref(x)-vWe(x, h, P), c(Tmin-10, Tmax), tol=0.001 )
    curve( vWe(x, h, P), xlim=c(pt[[1]], Tmax),  lwd=lw, add=TRUE )
    s <- ( vWe(Tmax-5, h, P)-vWe(Tmax, h, P) )/5                      # AV 5
    a <- 90-atan(s*sr)*180/pi
    lb <- format(h)
    text( x=pt[[1]], y=vWe(pt[[1]], h, P), label=lb, cex=0.6, srt=a )
}

if( VspNote ) {
    text( x=Tmin+Tstep.g, y=Wmax-3*Wstep.l, label="You can calculate specific volume using:", cex=0.6, pos=4 )
    text( x=Tmin+Tstep.g, y=Wmax-3.6*Wstep.l, label=expression(V[sp]==frac(0.287055*T,P)~(1+1.6078*W) ), cex=0.6, pos=4 )
}


dev.off()
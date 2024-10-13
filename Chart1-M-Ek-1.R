# Chart 1: M/Ek/1

# Obtaining ARL without randomization
inar1.arl <- function(L, U, rho, beta, k) {
  i    <- L:U
  d    <- U-L+1
  qij_ <- function(i, j) {
    if(i<=1) dnbinom(j, k, k/(k+rho)) else dnbinom(j-i+1, k, k/(k+rho))
  }
  qij  <- Vectorize(qij_)
  Q    <- outer(i, i, qij)
  one  <- rep(1, d)
  I    <- diag(1, d)
  cARL <- solve(I-Q, one)
  ARL  <- cARL[1]
  ARL
}


# Obtaining ARL with randomization
inar1.arl2 <- function(L, U, gL, gU, rho, beta, k) {
  i     <- L:U
  d     <- U-L+1
  qij_  <- function(i, j) {
    if(i<=1) dnbinom(j, k, k/(k+rho)) else dnbinom(j-i+1, k, k/(k+rho))
  }
  qij   <- Vectorize(qij_)
  Q     <- outer(i, i, qij)
  Q[,1] <- (1-gL) * Q[,1]
  Q[,d] <- (1-gU) * Q[,d]
  one   <- rep(1, d)
  I     <- diag(1, d)
  cARL  <- solve(I-Q, one)
  
  cARL[1]
}


# INAR(1) w/o randomization, L given, U searched for
inar1.get.U <- function(L, ARL0, rho, beta, k, OUTPUT=FALSE) {
  U   <- 1
  ARL <- inar1.arl(L, U, rho, beta, k)
  if ( OUTPUT ) cat(paste(U, "\t", ARL, "\t", ARL0, "\n"))
  dARL <- 1
  while ( ARL < ARL0 & dARL > 1e-6 ) {
    ARL.old <- ARL
    U   <- U + 1
    ARL <- inar1.arl(L, U, rho, beta, k)
    dARL <- abs(ARL - ARL.old)
    if ( OUTPUT ) cat(paste(U, "\t", ARL, "\t", ARL0, "\n"))
  }
  if ( ARL < ARL0 ) U <- NA
  U
}


# INAR(1) w/o randomization, (L,U) with minimal L so that ARL = ARL(rho) is increasing at given rho
inar1.get.UL <- function(ARL0, rho, beta, k, OUTPUT=FALSE) {
  L <- 0
  U <- inar1.get.U(L, ARL0, rho, beta, k, OUTPUT=OUTPUT)
  lARL <- inar1.arl(L, U, rho-1e-3, beta, k)
  ARL  <- inar1.arl(L, U, rho, beta, k)
  while ( ARL < lARL ) {
    L    <- L + 1
    U    <- inar1.get.U(L, ARL0, rho, beta, k, OUTPUT=OUTPUT)
    if ( is.na(U) ) {
      L    <- L - 1
      U    <- inar1.get.U(L, ARL0, rho, beta, k, OUTPUT=OUTPUT)
      break
    }
    lARL <- inar1.arl(L, U, rho-1e-2, beta, k)
    ARL  <- inar1.arl(L, U, rho, beta, k)
    if ( OUTPUT ) cat(paste(L, "\t", U, "\t", lARL, "\t", ARL, "\t\t", ARL0, "\n"))
  }
  c(L, U)
}


# INAR(1) w/ randomization, L, U, gL given, gU searched for
inar1.get.gU <- function(L, U, gL, ARL0, rho, beta, k, OUTPUT=FALSE) {
  minARL <- inar1.arl2(L, U, gL, 1, rho, beta, k)
  maxARL <- inar1.arl2(L, U, gL, 0, rho, beta, k)
  if ( OUTPUT ) cat(paste(minARL, "\t", maxARL, "\t\t", ARL0, "\n"))
  if ( minARL < ARL0 & ARL0 < maxARL ) {
    # starting values
    gU1  <- 1
    ARL1 <- minARL
    while ( ARL1 < ARL0 & gU1 > 0.1 ) {
      gU2  <- gU1
      ARL2 <- ARL1
      gU1  <- gU1 - .1
      ARL1 <- inar1.arl2(L, U, gL, gU1, rho, beta, k)
    }
    if ( gU1 < .05 ) {
      gU1  <- 0.05
      ARL1 <- inar1.arl2(L, U, gL, gU1, rho, beta, k)
    }
    # secant rule
    gU.error <- 1
    L.error  <- 1
    while ( gU.error > 1e-10  &  L.error > 1e-9 ) {
      gU3  <- gU1 + (ARL0 - ARL1)/(ARL2 - ARL1)*(gU2 - gU1)
      gU3 <- max(0, gU3)
      gU3 <- min(1, gU3)
      ARL3 <- inar1.arl2(L, U, gL, gU3, rho, beta, k)
      if ( OUTPUT ) cat(paste(gU3, "\t", ARL3, "\t\t", ARL0, "\n"))
      gU1 <- gU2; gU2 <- gU3
      ARL1 <- ARL2; ARL2 <- ARL3
      L.error  <- abs(ARL2 - ARL0)
      gU.error <- abs(gU2 - gU1)
    }
  } else {
    gU3 <- NA                      # old
    if ( minARL > ARL0 ) gU3 <- -2 # minimal reachable ARL too large
    if ( ARL0 > maxARL ) gU3 <- -1 # maximal reachable ARL too small
  }
  gU3
}


# INAR(1) w/ randomization, L, U, gU given, gL searched for
inar1.get.gL <- function(L, U, gU, ARL0, rho, beta, k, OUTPUT=FALSE) {
  minARL <- inar1.arl2(L, U, 1, gU, rho, beta, k)
  maxARL <- inar1.arl2(L, U, 0, gU, rho, beta, k)
  if ( OUTPUT ) cat(paste(minARL, "\t", maxARL, "\t\t", ARL0, "\n"))
  if ( minARL < ARL0 & ARL0 < maxARL ) {
    # starting values
    gL1  <- 1
    ARL1 <- minARL
    while ( ARL1 < ARL0 & gL1 > .1 ) {
      gL2  <- gL1
      ARL2 <- ARL1
      gL1  <- gL1 - .1
      ARL1 <- inar1.arl2(L, U, gL1, gU, rho, beta, k)
    }
    if ( gL1 < 0.05 ) {
      gL1  <- 0.05
      ARL1 <- inar1.arl2(L, U, gL1, gU, rho, beta, k)
    }
    # secant rule
    gL.error <- 1
    L.error  <- 1
    while ( gL.error > 1e-10  &  L.error > 1e-9 ) {
      gL3  <- gL1 + (ARL0 - ARL1)/(ARL2 - ARL1)*(gL2 - gL1)
      gL3 <- max(0, gL3)
      gL3 <- min(1, gL3)
      ARL3 <- inar1.arl2(L, U, gL3, gU, rho, beta, k)
      if ( OUTPUT ) cat(paste(gL3, "\t", ARL3, "\t\t", ARL0, "\n"))
      gL1 <- gL2; gL2 <- gL3
      ARL1 <- ARL2; ARL2 <- ARL3
      L.error  <- abs(ARL2 - ARL0)
      gL.error <- abs(gL2 - gL1)
    }
  } else {
    gL3 <- NA                      # old
    if ( minARL > ARL0 ) gL3 <- -2 # minimal reachable ARL too large
    if ( ARL0 > maxARL ) gL3 <- -1 # maximal reachable ARL too small
  }
  gL3
}


# some helpers
inar1.min.gL <- function(L, U, ARL0, rho, beta, k, OUTPUT=FALSE) {
  gL  <- 0
  gU  <- inar1.get.gU(L, U, gL, ARL0, rho, beta, k, OUTPUT=OUTPUT)
  if ( -1.5 < gU & gU < 0 ) gL <- NA
  if ( gU < -1.5 ) {
    for ( dig in 1:9 ) {
      if ( dig %% 2 == 1 ) {
        while ( gU < -1.5 & gL < 1-1e-10 ) {
          gL <- gL + 10^(-dig)
          gU <- inar1.get.gU(L, U, gL, ARL0, rho, beta, k, OUTPUT=OUTPUT)
          if ( OUTPUT ) cat(paste("gL =", gL, ",\tgU =", gU, "\n"))
        }
        if ( gU < -1.5 ) {
          gL <- NA
          break
        }
      } else {
        while ( -1.5 < gU  &  gL > 1e-10 ) {
          gL <- gL - 10^(-dig)
          gU <- inar1.get.gU(L, U, gL, ARL0, rho, beta, k, OUTPUT=OUTPUT)
          if ( OUTPUT ) cat(paste("gL =", gL, ",\tgU =", gU, "\n"))
        }
      }
    }
  }
  gL
}


inar1.min.gU <- function(L, U, ARL0, rho, beta, k, OUTPUT=FALSE) {
  gU  <- 0
  gL  <- inar1.get.gL(L, U, gU, ARL0, rho, beta, k, OUTPUT=OUTPUT)
  if ( -1.5 < gL & gL < 0 ) gU <- NA
  if ( gL < -1.5 ) {
    for ( dig in 1:9 ) {
      if ( dig %% 2 == 1 ) {
        while ( gL < -1.5 & gU < 1 ) {
          gU <- gU + 10^(-dig)
          gL <- inar1.get.gL(L, U, gU, ARL0, rho, beta, k, OUTPUT=FALSE)
          if ( OUTPUT ) cat(paste("gU =", gU, ",\tgL =", gL, "\n"))
        }
        if ( gL < -1.5 ) {
          gU <- NA
          break
        }
      } else {
        while ( -1.5 < gL & gU > 0 ) {
          gU <- gU - 10^(-dig)
          gL <- inar1.get.gL(L, U, gU, ARL0, rho, beta, k, OUTPUT=FALSE)
          if ( OUTPUT ) cat(paste("gU =", gU, ",\tgL =", gL, "\n"))
        }
      }
    }
  }
  gU
}



# INAR(1) w/ randomization, get all 4!
inar1.get.UL2 <- function(ARL0, rho, beta, k, target="rho", OUTPUT=FALSE, eps=1e-6, delta=1e-4, progress = NULL) {
  if (target == "rho") {
    r1 <- rho - eps
    r2 <- rho + eps
    b1 <- beta
    b2 <- beta
  }
  if (target == "beta") {
    r1 <- rho
    r2 <- rho
    b1 <- beta - eps
    b2 <- beta + eps
  }
  LU  <- inar1.get.UL(ARL0, rho, beta, k, OUTPUT=OUTPUT)
  L2  <- LU[1]
  U2  <- LU[2]
  if (L2 > 0) {
    L1  <- L2 - 1
    U1  <- inar1.get.U(L1, ARL0, rho, beta, k, OUTPUT=OUTPUT)
  } else {
    L1 <- 0
    U1 <- U2
    U2 <- U1 + 1
  }
  U2 <- U2 + 20
  
  # The time elapsed starts here
  total_steps <- (L2 - L1 + 1) * (U2 - U1 + 1)
  current_step <- 0
  
  start_time <- Sys.time()
  
  for (L in L1:L2) {
    for (U in U1:U2) {
      current_step <- current_step + 1
      elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      if (!is.null(progress)) {
        progress$inc(1 / total_steps, message = "Calculating limits", detail = paste("Elapsed time:", round(elapsed_time, 1), "seconds"))
      }
      
      gL1  <- inar1.min.gL(L, U, ARL0, rho, beta, k, OUTPUT=OUTPUT)
      if (is.na(gL1)) next
      gU1  <- inar1.get.gU(L, U, gL1, ARL0, rho, beta, k)
      lARL1 <- inar1.arl2(L, U, gL1, gU1, r1, b1, k)
      rARL1 <- inar1.arl2(L, U, gL1, gU1, r2, b2, k)
      dratio1 <- (rARL1 - lARL1) / (2 * eps)
      
      gU2 <- inar1.min.gU(L, U, ARL0, rho, beta, k, OUTPUT=OUTPUT)
      if (is.na(gU2)) next
      gL2 <- inar1.get.gL(L, U, gU2, ARL0, rho, beta, k)
      lARL2 <- inar1.arl2(L, U, gL2, gU2, r1, b1, k)
      rARL2 <- inar1.arl2(L, U, gL2, gU2, r2, b2, k)
      dratio2 <- (rARL2 - lARL2) / (2 * eps)
      
      if (OUTPUT) cat(paste("L =", L, ",\t U =", U, ",\tdr1 =", dratio1, ",\tdr2 =", dratio2, "\n"))
      
      if (dratio1 * dratio2 < 0) {
        # secant rule
        gL.error <- 1
        dr.error <- 1
        while (gL.error > 1e-10 & dr.error > delta) {
          gL3 <- gL1 + (0 - dratio1) / (dratio2 - dratio1) * (gL2 - gL1)
          gU3 <- inar1.get.gU(L, U, gL3, ARL0, rho, beta, k)
          lARL3 <- inar1.arl2(L, U, gL3, gU3, r1, b1, k)
          rARL3 <- inar1.arl2(L, U, gL3, gU3, r2, b2, k)
          dratio3 <- (rARL3 - lARL3) / (2 * eps)
          if (OUTPUT) cat(paste(gL3, "\t", dratio3, "\n"))
          gL1 <- gL2
          gL2 <- gL3
          dratio1 <- dratio2
          dratio2 <- dratio3
          dr.error  <- abs(dratio2)
          gL.error <- abs(gL2 - gL1)
        }
        L0 <- L
        U0 <- U
        L <- L2 + 1
        U <- U2 + 1
      }
      if (U > U2) break
    }
    if (L > L2) break
  }
  data.frame(L=L0, gL=gL3, U=U0, gU=gU3)
}


#########################################
#########################################
#########################################

# Load necessary libraries
library(shiny)
library(ggplot2)
library(shinycssloaders)

# Define UI ----
ui <- fluidPage(
  titlePanel("ARL-unbiased Xn-chart"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("rho",
                  HTML("Target Traffic Intensity (&rho;0):"),
                  min = 0.01,
                  max = 0.99,
                  value = 0.9,
                  step = 0.01),
      sliderInput("ARL0",
                  "Target ARL:",
                  min = 100,
                  max = 1000,
                  value = 200,
                  step = 1),
      sliderInput("k",
                  "Erlang Parameter (k):",
                  min = 1,
                  max = 100,
                  value = 5,
                  step = 1),
      downloadButton("downloadPlot", "Download Plot")
    ),
    
    mainPanel(
      withSpinner(plotOutput("ARLProfilePlot"), type = 8, color = "#0D6EFD")
    )
  )
)

# Define server logic ----
server <- function(input, output) {
  
  # Store the calculated values in a reactiveValues object to share between renderPlot and download
  rv <- reactiveValues()
  
  # Render the plot with progress bar and elapsed time
  output$ARLProfilePlot <- renderPlot({
    # Create a progress object
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    
    progress$set(message = "Calculating", value = 0)
    
    ARL0 <- input$ARL0
    rho0 <- input$rho
    k <- input$k
    
    # Increment the progress bar
    progress$inc(0.2, message = "Initializing variables", detail = "Elapsed time: calculating...")
    
    # Call the modified inar1.get.UL2 function with the progress object
    LUg <- inar1.get.UL2(ARL0, rho0, 0, k, OUTPUT = FALSE, progress = progress)
    
    # Store results in reactiveValues to use later in download
    rv$L3 <- LUg$L
    rv$gL <- LUg$gL
    rv$U3 <- LUg$U
    rv$gU <- LUg$gU
    
    progress$inc(0.3, message = "Generating plot data", detail = "Elapsed time: calculating...")
    
    rho <- seq(0.01, 0.99, by = 0.01)
    LL3 <- sapply(rho, function(r) inar1.arl2(rv$L3, rv$U3, rv$gL, rv$gU, r, 0, k))
    
    rv$rho <- rho  # Store rho and LL3 for download
    rv$LL3 <- LL3
    
    progress$inc(0.3, message = "Rendering plot", detail = "Elapsed time: calculating...")
    
    # Determine legend position based on rho0
    legend_pos <- if (rho0 > 0.5) "topleft" else "topright"
    
    par(mar = c(5, 5, 4, 2) + 0.1)
    
    # Plot
    plot(rho, LL3, type = 'l', col = 'black', xlab = "ρ", ylab = expression("ARL"))
    abline(v=rho0, h=ARL0, lty=4, col="grey")
    legend(legend_pos, legend = c(paste("LCL =", 0), paste("UCL =", round(rv$U3, 4)), 
                                  paste("γL =", round(rv$gL, 6)), paste("γU =", round(rv$gU, 6))), 
           cex = 0.8, inset = c(0.025, 0.1))
  })
  
  # Download the plot as an image
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("ARLProfilePlot", ".png", sep = "")
    },
    content = function(file) {
      png(file)
      
      # Reuse the cached plot data
      rho <- rv$rho
      LL3 <- rv$LL3
      U3 <- rv$U3
      gL <- rv$gL
      gU <- rv$gU
      ARL0 <- input$ARL0
      rho0 <- input$rho
      
      # Determine legend position based on rho0
      legend_pos <- if (rho0 > 0.5) "topleft" else "topright"
      
      par(mar = c(5, 5, 4, 2) + 0.1)
      
      # Plot for saving
      plot(rho, LL3, type = 'l', col = 'black', xlab = "ρ", ylab = expression("ARL"))
      abline(v=rho0, h=ARL0, lty=4, col="grey")
      legend(legend_pos, legend = c(paste("LCL =", 0), paste("UCL =", round(U3, 4)), 
                                    paste("γL =", round(gL, 6)), paste("γU =", round(gU, 6))), 
             cex = 0.8, inset = c(0.025, 0.1))
      
      dev.off()
    }
  )
}

# Run the app ----
shinyApp(ui = ui, server = server)

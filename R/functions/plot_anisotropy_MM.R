# Function to plot an anisotropy ellipse, given a fitted TMB model
# Original code from tinyVAST, modified by Markus Min

plot_anisotropy_MM <-
  function( Obj,
            FileName,
            nspecies,
            ControlList = list("Width"=4, "Height"=5, "Res"=200, "Units"='in'),
            type = "ellipse",
            Report = Obj$report(),
            TmbData = Obj$env$data ){
    
    # extract map
    Map = Obj$env$map
    Params = Obj$env$last.par.best
    
    # Decomposition
    Eigen = eigen(Report$H_b)
    
    # Arrows
    if( type=="arrow" ){
      png(file=FileName, width=ControlList$Width, height=ControlList$Height, res=ControlList$Res, units=ControlList$Units)
      par( mar=c(2,2,0,0), mgp=c(1.5,0.5,0), tck=-0.02)
      plot( 1, type="n", xlim=c(-1,1)*max(Eigen$values), ylim=c(-1,1)*max(Eigen$values))
      arrows( x0=rep(0,2), y0=rep(0,2), x1=Eigen$vectors[1,]*Eigen$values, y1=Eigen$vectors[2,]*Eigen$values)
      dev.off()
    }
    
    # Ellipses
    if( type=="ellipse" ){
      rss = function(V) sqrt(sum(V[1]^2+V[2]^2))
      Major_1 = Minor_1 = Major_2 = Minor_2 = NA
      # use estimated range to get these values
      Major_1 = Eigen$vectors[,1]*Eigen$values[1] * Report$Range_b
      Minor_1 = Eigen$vectors[,2]*Eigen$values[2] * Report$Range_b
      
      png(file=FileName, width=ControlList$Width, height=ControlList$Height, res=ControlList$Res, units=ControlList$Units)
      par( mar=c(3,3,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
      Range = 1.1 * c(-1,1) * max(abs( cbind(Major_1,Minor_1, Major_2,Minor_2) ),na.rm=TRUE)
      plot( 1, type="n", xlim=Range, ylim=c(Range[1],Range[2]*1.2), xlab="", ylab="")
      # plot the ellipse
      shape::plotellipse( rx=rss(Major_1), ry=rss(Minor_1), angle=-1*(atan(Major_1[1]/Major_1[2])/(2*pi)*360-90), lcol="green", lty="solid")
      title( "Distance at 10% correlation" )
      mtext(side=1, outer=FALSE, line=2, text="Eastings (km.)")
      mtext(side=2, outer=FALSE, line=2, text="Northings (km.)")
      legend( "top", legend=c("1st linear predictor"), fill=c("green"), bty="n")
      #abline( h=0, v=0, lty="dotted")
      dev.off()
    }
  }
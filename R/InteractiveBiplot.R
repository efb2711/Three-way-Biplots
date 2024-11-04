#' Interactive Biplot
#'
#' Generates an interactive biplot based on a three-dimensional array `Y3D` and
#' its Tucker3 or PARAFAC decomposition.
#'
#' @param Y3D The unfolding of the first mode of the three-way data array.
#' @param G The unfolding of the first mode of the core matrix (only for
#' Tucker3 decomposition)
#' @param A Component matrix of the first mode
#' @param B Component matrix of the second mode
#' @param C Component matrix of the third mode
#' @param namesA Optional vector of row names for `A`. If not provided,
#' default names are generated.
#' @param namesBC Optional vector of column names for the Kronecker product of
#' `B` and `C`. If not provided, default names are generated based on `J` and
#' `K` dimensions.
#'
#' @return A list with class `ContinuousBiplot`
#'
#' @examples
#' # Example of function usage:
#' data(X3way)
#' Y3D <- X3way
#' data(ResultsTucker3)
#' namesA = c("i1","i2","i3","i4","i5","i6","i7","i8")
#' Biplot <- InteractiveBiplot(Y3D = Y3D, A = A, B = B, C = C, G = G,
#' namesA = namesA)
#' colores <- c("red", "blue", "dark green", "purple", "orange", "pink","brown")
#' vector_colores <- rep(colores, 7)
#' plot(Biplot,ColorVar = vector_colores,mode="ah",ColorInd = "black")
#'
#' @export
#'
InteractiveBiplot <- function(Y3D, A, B, C, G = NULL,
                              namesA=NULL, namesBC=NULL)
{
  r1 = ncol(A)
  I = nrow(A)
  r2 = ncol(B)
  J = nrow(B)
  r3 = ncol(C)
  K = nrow(C)
  if (is.null(namesA)) namesA= paste("A",1:I)
  if (is.null(namesBC)) {
    names=NULL
    for (k in 1:K)
      names=c(names, paste("S",1:J,"-R",k, sep=""))
    colnames(Y3D)=names
    namesBC = names
  }

  Biplot=list()
  Biplot$Title = "Interactive Biplot"
  Biplot$Type = "Interactive"
  Biplot$Non_Scaled_Data = Y3D
  Biplot$alpha=0
  Biplot$Dimension=min(r1, r2, r3)
  Biplot$Means = apply(Y3D, 2, mean)
  Biplot$Medians = apply(Y3D, 2, median)
  Biplot$Deviations = apply(Y3D, 2, sd)
  Biplot$Minima = apply(Y3D, 2, min)
  Biplot$Maxima = apply(Y3D, 2, max)
  Biplot$P25 = apply(Y3D, 2, quantile)[2, ]
  Biplot$P75 = apply(Y3D, 2, quantile)[4, ]
  Biplot$Scaled_Data = Y3D

  Biplot$nrows = I
  Biplot$ncols = J*K
  Biplot$dim = r1

  sct=sum(Y3D^2)
  scf = apply((Y3D^2), 1, sum)
  scc = apply((Y3D^2), 2, sum)

  rownames(A) = namesA
  if (is.null(G)) {
    BC = khatri_rao( C , B)
  }else {
    BC = t(G %*% kronecker( t(C) , t(B) ))
  }
  Esp=A%*%t(BC)

  sce=rep(0,r1)
  cf=matrix(0,I,r1)
  cc=matrix(0,J*K,r1)
  for (i in 1:r1){
    Yesp=matrix(A[,i])%*%t(matrix(BC[,i]))
    scfe = apply((Yesp^2), 1, sum)
    scce = apply((Yesp^2), 2, sum)
    cf[,i]=scfe/scf
    cc[,i]=scce/scc
    sce[i]=sum(Yesp^2)
  }

  cf=100*cf
  cc=100*cc
  Biplot$EigenValues = sce
  Biplot$Inertia = 100 * sce/sct

  Biplot$CumInertia = cumsum(Biplot$Inertia)
  sca = sum(A^2)
  scb = sum(BC^2)
  sca = sca/I
  scb = scb/(J*K)
  scf = sqrt(sqrt(scb/sca))
  A = A * scf
  BC = BC/scf

  rownames(BC) = namesBC
  colnames(A) = paste("Dim", 1:r1)

  Biplot$RowCoordinates=A
  Biplot$ColCoordinates=BC


  rownames(cf) = namesA
  colnames(cf) = paste("Dim", 1:r1)
  Biplot$RowContributions=cf

  rownames(cc) = namesBC
  colnames(cc) = paste("Dim", 1:r1)
  Biplot$ColContributions=cc

  Biplot$Scale_Factor = scf

  class(Biplot)="ContinuousBiplot"
  return (Biplot)
}


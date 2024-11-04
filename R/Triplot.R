#' Triplot
#'
#' Generates a triplot based on a three-dimensional array `Y3D`, a data
#' matrix X and their decomposition.
#'
#' @param Y3D The unfolding of the first mode of the three-way data array.
#' @param X A data matrix
#' @param G The unfolding of the first mode of the core matrix
#' @param A Component matrix of the first mode
#' @param B1 Component matrix of the second mode
#' @param C Component matrix of the third mode
#' @param B2 Covariates matrix
#' @param namesA Optional vector of row names for `A`. If not provided,
#' default names are generated.
#' @param namesBC Optional vector of column names for the Kronecker product of
#' `B` and `C`. If not provided, default names are generated based on `J` and
#' `K` dimensions.
#' @param namesB2 Optional vector of row names for `B2`. If not provided,
#' default names are generated.
#' @return A list with class `ContinuousBiplot`
#'
#' @examples
#' # Example of function usage:
#' data(X3way)
#' Y3D <- X3way
#' data(X2way)
#' X <- X2way
#' data(ResultsTucker3)
#' namesA = c("i1","i2","i3","i4","i5","i6","i7","i8")
#' namesB2 = c("Fear to be refused","Kindness",
#' "Importance of others’ judgments","Altruism","Neuroticism","Openness",
#' "Being strict to oneself","Low selfesteem","Conscientiousness","Depression")
#' colores <- c("red", "blue", "dark green", "purple", "orange", "pink",
#' "brown")
#' vector_colores <- rep(colores, 7)
#' Triplot <- Triplot (X = X,Y3D=Y3D, A = A, B1 = B1, C = C, B2 = B2,
#'                    namesA = namesA, namesB2 = namesB2)
#' plot(Triplot,ColorVar = vector_colores,mode="ah",ColorInd = "black")
#'
#' @export
#'
Triplot <- function(X,Y3D, G = NULL, A, B1, C, B2,
                              namesA=NULL, namesBC=NULL,namesB2=NULL)
{
  # Pinta el biplot interactivo
  BiplotY <- InteractiveBiplot(Y3D=Y3D,A=A,B=B1,C=C,G=G,namesA = namesA,
                               namesBC=namesBC)

  # Construcción del biplot para X (separado)
  I = nrow(A)
  L = nrow(B2)
  r = ncol(A)

  BiplotX = list()
  BiplotX$nrows = I
  BiplotX$ncol = L
  BiplotX$dim = r

  sct=sum(X^2)

  scf = apply((X^2), 1, sum)
  scc = apply((X^2), 2, sum)

  sce=rep(0,r)
  cf=matrix(0,I,r)
  cc=matrix(0,L,r)
  for (i in 1:r){
    Xesp=matrix(A[,i])%*%t(matrix(B2[,i]))
    scfe = apply((Xesp^2), 1, sum)
    scce = apply((Xesp^2), 2, sum)
    cf[,i]=scfe/scf
    cc[,i]=scce/scc
    sce[i]=sum(Xesp^2)
  }
  cf=100*cf
  cc=100*cc
  BiplotX$EigenValues = sce
  BiplotX$Inertia = 100 * sce/sct

  BiplotX$CumInertia = cumsum(BiplotX$Inertia)
  #BiplotX$EV = t3pc$B2

  sca = sum(A^2)
  scb = sum(B2^2)
  sca = sca/I
  scb = scb/L
  scf = sqrt(sqrt(scb/sca))
  a = A * scf
  b = B2/scf

  rownames(a) <- namesA
  rownames(b) <- namesB2
  BiplotX$RowCoordinates=a
  BiplotX$ColCoordinates=b



  rownames(cf) = namesA
  colnames(cf) = paste("Dim", 1:r)
  BiplotX$RowContributions=cf

  rownames(cc) = colnames(X)
  colnames(cc) = paste("Dim", 1:r)
  BiplotX$ColContributions=cc

  class(BiplotX)="ContinuousBiplot"

  BiplotY$ContSupVarsBiplot=BiplotX
  class(BiplotY$ContSupVarsBiplot)="ContSupVarsBiplot"

  return(BiplotY)
}

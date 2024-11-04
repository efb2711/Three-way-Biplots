#' Joint Biplot
#'
#' Generates a joint biplot based on a three-dimensional array `Y3D` and
#' its Tucker3 or PARAFAC decomposition.
#'
#' @param A Component matrix of the first mode
#' @param B Component matrix of the second mode
#' @param C Component matrix of the third mode
#' @param G The unfolding of the first mode of the core matrix (only for
#' Tucker3 decomposition)
#' @param fixmode Mode to be set. Possible values 1, 2 or 3
#' @param fixunit Component number to be setComponent number to be set
#' @param namesA Optional vector of row names for `A`. If not provided,
#' default names are generated.
#' @param namesB Optional vector of row names for `B`. If not provided,
#' default names are generated.
#' @param namesC Optional vector of row names for `C`. If not provided,
#' default names are generated.
#'
#' @return A list with class `ContinuousBiplot`
#'
#' @examples
#' # Example of function usage:
#' data(ResultsTucker3)
#' namesA = c("i1","i2","i3","i4","i5","i6","i7","i8")
#' namesB = c("Other anger","Shame","Love","Sorrow","Fear","Guilt",
#' "Self(anger)")#'
#' namesC = c("Quarrelling with someone","Partner leaves you",
#' "Someone is telling lies about you","Giving a bad speech","Failing a test",
#' "Writing a bad paper")
#' Biplot <- JointBiplot(A = A, B = B, C = C, G = G, fixmode = 3, fixunit = 1,
#' namesA = namesA, namesB = namesB, namesC = namesC)
#' colores <- c("red", "blue", "dark green", "purple", "orange", "pink",
#' "brown")
#' vector_colores <- rep(colores, 7)
#' plot(Biplot,ColorVar = vector_colores,mode="ah",ColorInd = "black")
#'
#' @export
JointBiplot <- function(A, B, C, G = NULL, fixmode, fixunit,
                         namesA=NULL, namesB=NULL, namesC=NULL)
{
  r1 = ncol(A)
  I = nrow(A)
  r2 = ncol(B)
  J = nrow(B)
  r3 = ncol(C)
  K = nrow(C)
  if (is.null(namesA)) namesA= paste("A",1:I)
  if (is.null(namesB)) namesB= paste("B",1:J)
  if (is.null(namesC)) namesC= paste("C",1:K)

  if (fixmode == 1) {
    if (is.null(G)) {
       Bmat = B
       Cmat = C
    } else {
      Gmat = matrix(t(G[fixunit, ]), r2, r3)
      Smat = B %*% Gmat %*% t(C)
      SVD = svd(Smat)
      Bmat = (J/K)^0.25 * SVD$u[,1:r2] %*% diag(SVD$d[1:r2])^0.5
      Cmat = (K/J)^0.25 * SVD$v[,1:r3] %*% diag(SVD$d[1:r3])^0.5
    }
    colnames(Bmat)=paste("Dim", 1:r2)
    colnames(Cmat)=paste("Dim", 1:r3)
    rownames(Bmat)=namesB
    rownames(Cmat)=namesC

    Biplot=list()
    Biplot$Type="JointBiplot"
    Biplot$Inertia = c(50,50)
    Biplot$CumInertia = cumsum(Biplot$Inertia)
    Biplot$RowCoordinates=Bmat
    Biplot$ColCoordinates=Cmat

    Biplot$RowContributions=matrix(0.5, J, 2)
    Biplot$ColContributions=matrix(0.5, K, 2)

    class(Biplot)="ContinuousBiplot"
  }
  if (fixmode == 2){
    if (is.null(G)) {
      Amat = A
      Cmat = C
    } else {
      GG = permnew(K, r1, r2, r3)
      Gmat = matrix(t(GG[fixunit, ]), r3, r1)
      Smat = C %*% Gmat %*% t(A)
      SVD = svd(Smat)
      Cmat = (K/I)^0.25 * SVD$u[, 1:r3] %*% diag(SVD$d[1:r3])^0.5
      Amat = (I/K)^0.25 * SVD$v[, 1:r1] %*% diag(SVD$d[1:r1])^0.5
    }

    rownames(Amat)=namesA
    rownames(Cmat)=namesC

    Biplot=list()
    Biplot$Type="JointBiplot"
    Biplot$alpha = 0.5
    Biplot$Inertia = c(50,50)
    Biplot$CumInertia = cumsum(Biplot$Inertia)
    Biplot$RowCoordinates=Amat
    Biplot$ColCoordinates=Cmat
    Biplot$RowContributions=matrix(0.5, I, 2)
    Biplot$ColContributions=matrix(0.5, K, 2)
    class(Biplot)="ContinuousBiplot"
  }
  if (fixmode == 3){
    if (is.null(G)) {
      Amat = A
      Bmat = B
    } else {
      GG = permnew(G, r1, r2, r3)
      GG = permnew(G, r2, r3, r1)
      Gmat = matrix(t(GG[fixunit, ]), r1, r2)
      Smat = A %*% Gmat %*% t(B)
      SVD = svd(Smat)

      Amat = (J/I)^0.25 * SVD$u[,1:r1] %*% diag(SVD$d[1:r1])^0.5
      Bmat = (I/J)^0.25 * SVD$v[,1:r2] %*% diag(SVD$d[1:r2])^0.5
    }
    sca = sum(Amat^2)
    scb = sum(Bmat^2)
    sca = sca/I
    scb = scb/J
    scf = sqrt(sqrt(scb/sca))
    Amat = Amat * scf
    Bmat = Bmat / scf

    colnames(Amat)=paste("Dim", 1:r1)
    colnames(Bmat)=paste("Dim", 1:r2)
    rownames(Amat)=namesA
    rownames(Bmat)=namesB

    Biplot=list()
    Biplot$Type="JointBiplot"
    Biplot$alpha = 0.5
    Biplot$Inertia = c(50,50)
    Biplot$CumInertia = cumsum(Biplot$Inertia)
    Biplot$RowCoordinates=Amat
    Biplot$ColCoordinates=Bmat
    Biplot$RowContributions=matrix(50, I, 2)
    Biplot$ColContributions=matrix(50, J, 2)
    class(Biplot)="ContinuousBiplot"
  }
  return(Biplot)
}


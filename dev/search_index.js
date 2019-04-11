var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "PosDefManifold Documentation",
    "title": "PosDefManifold Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "#PosDefManifold-Documentation-1",
    "page": "PosDefManifold Documentation",
    "title": "PosDefManifold Documentation",
    "category": "section",
    "text": ""
},

{
    "location": "#Installation-1",
    "page": "PosDefManifold Documentation",
    "title": "Installation",
    "category": "section",
    "text": ""
},

{
    "location": "#Overview-1",
    "page": "PosDefManifold Documentation",
    "title": "Overview",
    "category": "section",
    "text": "(Image: Figure 1)"
},

{
    "location": "#Code-units-1",
    "page": "PosDefManifold Documentation",
    "title": "Code units",
    "category": "section",
    "text": "PosDefManifold includes five code units (.jl files):Unit Description\nPosDefManifold.jl Main module\nriemannianGeometry.jl The fundamental unit collecting all Riemannian functions\nlinearAlgebra.jl Collection of linear algebra routines\nsignalProcessing.jl Collection of signal processing routines\ntest.jl Unit performing all tests"
},

{
    "location": "#Contents-1",
    "page": "PosDefManifold Documentation",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [ \"index.md\", \"introToRiemannianGeometry.md\", \"MainModule.md\", \"riemannianGeometry.md\", \"linearAlgebra.md\", \"signalProcessing.md\", \"test.md\"]\nDepth = 1"
},

{
    "location": "#Index-1",
    "page": "PosDefManifold Documentation",
    "title": "Index",
    "category": "section",
    "text": ""
},

{
    "location": "introToRiemannianGeometry/#",
    "page": "Intro to Riemannian Geometry",
    "title": "Intro to Riemannian Geometry",
    "category": "page",
    "text": ""
},

{
    "location": "introToRiemannianGeometry/#Intro-to-Riemannian-Geometry-1",
    "page": "Intro to Riemannian Geometry",
    "title": "Intro to Riemannian Geometry",
    "category": "section",
    "text": "The study of appropriate distance measures for positive definite matrices has recently grown very fast, driven by practical problems in radar data processing, image processing, computer vision, shape analysis, medical imaging (especially diffusion MRI and Brain-Computer Interface), sensor networks, elasticity, mechanics, numerical analysis and machine learning (e.g., see references in Congedo et al., 2017a)🎓.In many applications the observed data can be conveniently summarized by positive definite matrices, which are either symmetric positive definite (SPD: real) or Hermitian Positive Definite (HPD: complex). For example, those may be some form of the data covariance matrix in the time, frequency or time-frequency domain, or autocorrelation matrices, kernels, slices of tensors, density matrices, elements of a search space, etc. Positive definite matrices are naturally treated as points on a smooth Riemannian manifold allowing useful operations such as interpolation, smoothing, filtering, approximation, averaging, signal detection and classification. Such operations are the object of the present PosDefManifold library.More formally, this Julia library treats operations on the metric space (P δ^2) of n・n positive definite matrices endowed with a distance or symmetric divergence delta(P x P)0 . Several matrix distances or matrix divergences δ are considered. Using some of them, the most important one being the Fisher metric, we define a Riemannian manifold. In mathematics, this is the subject of Riemannian geometry and information geometry.Note that throughout this library the word \'metric\' is used loosely for referring to the actual Riemannian metric on the tangent space and to the resulting distance or to general symmetric divergence acting on P, regardless the fact that we are dealing with a metric in the strict sense and that it induces or not a Riemannian geometry in P. This is done for convenience of exposition, since in practice those \'metrics\' in PosDefManifold may be used interchangeably."
},

{
    "location": "introToRiemannianGeometry/#Riemannian-manifolds-1",
    "page": "Intro to Riemannian Geometry",
    "title": "Riemannian manifolds",
    "category": "section",
    "text": "Here are some important definitions:A smooth manifold in differential geometry is a topological space that is locally similar to the Euclidean space and has a globally defined differential structure.The tangent space at point G is the vector space containing the tangent vectors to all curves on the manifold passing through G (Fig. 1).A smooth Riemannian manifold is equipped with an inner product on the tangent space (a Riemannian metric) defined at each point and varying smoothly from point to point. For manifold P the tangent space is the space of symmetric or Hermitian matrices.Thus, a Riemannian metric turns the metric space (P δ^2) into a Riemannian manifold. This is the case, for example, of the Fisher metric, which has a fundamental role in the manifolds of positive definite matrices and of the Wasserstein metric, fundamental in optimal transport theory.(Image: ) Figure 1"
},

{
    "location": "introToRiemannianGeometry/#geodesic-1",
    "page": "Intro to Riemannian Geometry",
    "title": "geodesic",
    "category": "section",
    "text": "The key object in the P manifold is the geodesic, the shortest path joining two points P and Q on the manifold, analogous to straight lines in the Euclidean space (Fig. 1). The gedesic equation with arclength 0a1 is the equation of the points along the path, where with a=0 we stay at P and with a=1 we move all the way to Q. The points along the geodesic in between P and Q (0a1) can be understood as weighted means of P and Q. For example, the geodesic equation according to the Euclidean metric is (1-a)P + aQ, which is the traditional way to define weighted means. With the metrics we consider here, geodesics are unique and always exist. Furthermore, as we will see, using the Fisher metric those geodesics extends indefinitely, i.e., they are definied and always remain positive definite for -a."
},

{
    "location": "introToRiemannianGeometry/#distance-1",
    "page": "Intro to Riemannian Geometry",
    "title": "distance",
    "category": "section",
    "text": "The length of the geodesic (at constant velocity) between two points gives the distance δ(P Q)."
},

{
    "location": "introToRiemannianGeometry/#distance-from-the-origin-1",
    "page": "Intro to Riemannian Geometry",
    "title": "distance from the origin",
    "category": "section",
    "text": "In contrast to an Euclidean space, the origin of the P manifold endowed with the Fisher metric is not 0_n, but I_n, the identity matrix of dimension n・n. The distance between a point P and the origin, i.e., δ(P I), is analogous therein to the length of vectors in Euclidean space. This Riemannian manifold is symmetric around I_n, i.e., δ(P I)=δ(P^-1 I) and δ(P Q)=δ(P^-1 Q^-1). This will be made more precise when we talk about invariances."
},

{
    "location": "introToRiemannianGeometry/#mean-1",
    "page": "Intro to Riemannian Geometry",
    "title": "mean",
    "category": "section",
    "text": "The mid-point on the geodesic relying P and Q is named the mean. Using the Euclidean metric this is the arithmetic mean of P and Q and using the inverse Euclidean metric this is their harmonic mean. As we will see, those are straightforward extensions of their scalar counterparts. Using the Fisher metric the mid-point of the geodesic relying P and Q allows the proper generalization to matrices of the scalars\' geometric mean. The other metrics allows other definition of means (see below)."
},

{
    "location": "introToRiemannianGeometry/#Fréchet-mean-1",
    "page": "Intro to Riemannian Geometry",
    "title": "Fréchet mean",
    "category": "section",
    "text": "Using Fréchet\'s variational approach we can extend to positive-definite matrices the concept of weighted mean of a set of scalars; as the midpoint G on the geodesic relying P and Q is the minimizer of sigma^2(P G)+sigma^2(Q G), so the mean G of points P_1 P_2P_k is the matrix G verifyingtextrmargmin_Gsum_i=1^kδ^2(P_iG)Thus, every metric induces a distance (or divergence) function, which, in turn, induces a mean."
},

{
    "location": "introToRiemannianGeometry/#invariances-1",
    "page": "Intro to Riemannian Geometry",
    "title": "invariances",
    "category": "section",
    "text": "An important characteristic of metrics is that they may induce invariance properties on the distance, which are in turn inherited by the mean.Let us denote shortly by P_i the set P_1P_k, where i=1k  and by GP_i the Fréchet mean of the set (in this section we drop the weights here for keeping the notation short). The most important invariance properties are:invariance effect on distance δ(PQ) effect on mean GP_i\nrotation δ(PQ)=δ(U^*PUU^*QU) GU^*P_iU=U^*GP_iU\naffinty δ(PQ)=δ(B^*PB^*B^*QB) GB^*P_iB=B^*GP_iB\ninversion δ(PQ)=δ(P^-1Q^-1) GP_i^-1=G^-1P_ifor any unitary U unitary and non-singular B.The affine invariance implies the rotation invariance and is also named congruence invariance."
},

{
    "location": "introToRiemannianGeometry/#metrics-1",
    "page": "Intro to Riemannian Geometry",
    "title": "metrics",
    "category": "section",
    "text": "We are interested in distance or divergence functions, the difference between the two being that a divergence does not need to be symmetric nor to satisfy the triangle inequality. Note that in PosDefManifold we consider only distances and symmetric divergences. In fact those are of greater interest in practice. One can find several distances and divergences in the literature and they often turn out to be related to each other, see for example (Chebby and Moakher, 2012; Cichocki et al., 2015; Sra, 2016)🎓. Ten of them are implemented in PosDefManifold and two of them are Riemannian metrics (the Fisher and Wasserstein metric as we have said). In this section we give a complete list of the expressions for their induceddistance of a point P from the origin,  \ndistance between two points P and Q,\ngeodesic relying P to Q (hence the weighted means of P and Q)\nweighted Fréchet mean G(Pw) of a set of k2 points P_1P_k with associated real non-negative weights w_1w_k summing up to 1.note: Nota Bene\nIn the following, the weights w_1w_k are always supposed summing up to 1, superscript * indicate conjugate transpose (or just transpose if the matrix is real) and if a is the arclength of a geodesic, we define for convenience b=1-a."
},

{
    "location": "introToRiemannianGeometry/#Euclidean-1",
    "page": "Intro to Riemannian Geometry",
    "title": "Euclidean",
    "category": "section",
    "text": "This is the classical Euclidean distance leading to the usual arithmetic mean. In general this metric is not well adapted to the P manifold. It verifies only the rotation invariance, however the mean also verifies the congruence invariance.distance² to I distance²\nP-I^2 P-Q^2geodesic Fréchet mean\nbP + aQ sum_i=1^kw_i P_i"
},

{
    "location": "introToRiemannianGeometry/#inverse-Euclidean-1",
    "page": "Intro to Riemannian Geometry",
    "title": "inverse Euclidean",
    "category": "section",
    "text": "This is the classical harmonic distance leading to the harmonic mean. It verifies only the rotation invariance, however the mean also verifies the congruence invariance.distance² to I distance²\nP^-1-I^2 P^-1-Q^-1^2geodesic Fréchet mean\nbig(bP^-1 + aQ^-1big)^-1 big(sum_i=1^kw_i P_i^-1big)^-1"
},

{
    "location": "introToRiemannianGeometry/#Cholesky-Euclidean-1",
    "page": "Intro to Riemannian Geometry",
    "title": "Cholesky Euclidean",
    "category": "section",
    "text": "This is a very simple metric that has been tried to improve the Euclidean one. It is rarely used (see for example Dai et al., 2016)🎓. It does not verify any invariance. Let L_P be the lower triangular Cholesky factor of P, thendistance² to I distance²\nL_P-I^2  L_P-L_Q ^2geodesic Fréchet mean\n(bL_P+aL_Q)(bL_P+aL_Q)^* big(sum_i=1^kw_i L_P_ibig)big(sum_i=1^kw_i L_P_ibig)^*"
},

{
    "location": "introToRiemannianGeometry/#log-Euclidean-1",
    "page": "Intro to Riemannian Geometry",
    "title": "log Euclidean",
    "category": "section",
    "text": "If matrices P_1P_k all pair-wise commute, then this metric coincides with the Fisher metric. See (Arsigny et al., 2007 ; Bhatia et al., 2019a)🎓. It enjoys the rotation and inversion invariance. The log-Euclidean distance to I is the same as per the Fisher metric. This mean has the same determinant as the Fisher mean, and trace equal or superior to the trace of the Fisher mean. A minimum trace log Euclidean mean approximating well the Fisher mean has been proposed in Congedo et al. (2015)🎓.distance² to I distance²\ntextrmlog(P)^2 textrmlog(P)-textrmlog(Q)^2geodesic Fréchet mean\ntextrmexpbig(textrmlogP + atextrmlogQbig) textrmexpbig(sum_i=1^kw_ihspace1pttextrmlogP_ibig)"
},

{
    "location": "introToRiemannianGeometry/#log-Cholesky-1",
    "page": "Intro to Riemannian Geometry",
    "title": "log Cholesky",
    "category": "section",
    "text": "It is a recently proposed distance in P. Like the Cholesky Euclidean metric here above, it exploits the diffeomorphism between matrices in P and their Cholesky factor, such that L_PL_P^*=P, thanks to the fact that the Cholesky factor is unique and that the map is smooth (Lin, 2019)🎓. The mean has the same determinant as the Fisher and log-Euclidean mean.Let L_X,S_X and D_X be the lower triangle, the strictly lower triangle and the diagonal part of X, respectively (hence, S_X+D_X=L_X), thenDistance² to I Distance²\nS_P-I^2+textrmlogD_P^2 S_P-S_Q^2+textrmlogD_P-textrmlogD_Q^2geodesic: S_P+a(S_Q-S_P)+D_Phspace2pttextrmexpbig(atextrmlogD_Q-atextrmlogD_Pbig)Fréchet mean: TT^*, where T=sum_i=1^kw_iS_P_i+sum_i=1^kw_itextrmlogD_P_i"
},

{
    "location": "introToRiemannianGeometry/#Fisher-1",
    "page": "Intro to Riemannian Geometry",
    "title": "Fisher",
    "category": "section",
    "text": "The Fisher metric, also known as affine-invariant, natural and Fisher-Rao metric, among others names, has a paramount importance for the P manifold, standing out as the natural choice both from the perspective of differential geometry and information geometry. Endowed with the Fisher metric the manifold P is Riemannian, has nonpositive curvature and is symmetric. This metric verifies all three invariances we have considered.Distance² to I Distance²\ntextrmlog(P)^2 textrmlog(P^-12QP^-12)^2geodesic\nP^12 big(P^-12 Q P^-12big)^a P^12Fréchet mean: it does not have a closed-form solution in general. The solution is the unique positive definite matrix G satisfying (Bhatia and Holbrook, 2006; Moakher, 2005).🎓sum_i=1^kw_itextrmlogbig(G^-12 P_i G^-12big)=0For estimating it, PosDefManifold implements a dedicated version of the MPM algorithm (Congedo et al., 2017b)🎓.This mean is known under many different names (Fisher, Rao, Fisher-Rao, Pusz-Woronowicz, Cartan, Fréchet, Karcher, geometric....).   The ‘centrality’ of this mean among a wide family of divergence-based means can be appreciated in Fig. 4 of Cichocki et al. (2015)🎓.The geometric mean G of two matrices P and Q is denoted shortly as PtextrmQ. Currently it is an object of intense study because of its interesting mathematical properties. For instance,it is the unique solution to Riccati equation GQ^-1G=P\nit is equal to F^-*D_1^12D_2^12F^-1 for whatever joint diagonalizer F of P and Q, i.e., for whatever matrix F satisfying F^*PF=D_1 and F^*QF=D_2, with D_1, D_1 non-singular diagonal matrices (Congedo et al., 2015)🎓.\nit enjoys all 10 properties of means postulated in the seminal work of Ando et al. (2010)🎓.When P and Q commutes, the Fisher mean of two matrices reduces to P^12Q^12, which indeed in this case is the log-Euclidean mean frac12textrmlogP + frac12textrmlogQ.The Fisher geodesic equation is usually denoted Ptextrm_aQ. Note that Itextrm_aP=P^a and Ptextrm_aI=P^b, where b=1-a.Fisher geodesic equation verifies Ptextrm_aQ=Qtextrm_bP and (Ptextrm_aQ)^-1=P^-1textrm_aQ^-1.An interesting property of the Fisher metric is that using its geodesic equation we can extrapolate positive matrices, always remaining in P. That is, using any real value of a :with 0  a  1 we move toward Q		(attraction),\nwith a  1 we move over and beyond Q	(extrapolation) and\nwith a 0 we move back away from Q 	(repulsion).Something similar can be done using the log Cholesky metric as well."
},

{
    "location": "introToRiemannianGeometry/#power-means-1",
    "page": "Intro to Riemannian Geometry",
    "title": "power means",
    "category": "section",
    "text": "The arithmetic, harmonic and geometric mean we have encountered are all members of the 1-parameter family of power means (with parameter p-1 1) introduced by Lim and Palfia (2012)🎓 to generalize the concept of power means of scalars (also known as Hölder means or generalized means). The family of power means G with parameter p satisfies equationG=sum_i=1^kw_ibig(Gtextrm_pP_ibig),whereGtextrm_pP_i is the Fisher geodesic equation we have discussed here above talking about the Fisher metric. In particular:with p=-1 this is the harmonic mean (see the inverse Euclidean metric)\nwith p=+1 this is the arithmetic mean (see the Euclidean metric)\nat the limit of p evaluated at zero from both side this is the geometric mean (see the Fisher metric).Thus, the family of power means continuously interpolate between the arithmetic and harmonic mean passing through the the geometric mean.All power means enjoy the congruence invariance (hence the rotation invariance), but only the geometric mean enjoy also the inversion invariance.The power mean with p=frac12 is the solution of the Fréchet mean problem using the following divergence (Bhatia, Gaubert and Jain, 2019)🎓δ^2(PQ)=textrmtr(P+Q)-2textrmtrPtextrmQ = textrmtr(textrmarithm mean(P Q))  textrmtr(textrmgeom mean(P Q))"
},

{
    "location": "introToRiemannianGeometry/#generalized-means-1",
    "page": "Intro to Riemannian Geometry",
    "title": "generalized means",
    "category": "section",
    "text": "When the matrices in the set all pairwise commute, it has been proved in Lim and Palfia (2012, see Property 1, p. 1502) 🎓 that the power means we have just seen reduce tobig(sum_i=1^kw_iP_i^pbig)^1p,which are the straightforward extension of scalar power means (see generalized means) to matrices. As usual, such straightforward extensions work well in commuting algebra, but not in general. See for example the case of the mean obtained using the log Euclidean metric, which is the straightforward extension to matrices of the scalar geometric mean, but is not the matrix geometric mean, unless the matrices all pairwise commute.Both the generalized means and the power means have a parameter p-1 1. For the latter, the solution is implemented via the fixed-point MPM algorithm (Congedo et al., 2017b)🎓."
},

{
    "location": "introToRiemannianGeometry/#modified-Bhattacharyya-mean-1",
    "page": "Intro to Riemannian Geometry",
    "title": "modified Bhattacharyya mean",
    "category": "section",
    "text": "If matrices P_1 P_2P_k all pair-wise commute, the special case p=frac12 yields the following instance of power means:big(sum_i=1^kw_iP_i^12big)^12.This mean has been proposed  in a different context by Moakher (2012)🎓 as a modified Bhattacharyya mean, since it is a modification of the Bhattacharyya mean we will encounter next under the name logdet zero. It is worth noting that in commuting algebra Moakher’s mean also corresponds to the mean obtained with the Wasserstein metric."
},

{
    "location": "introToRiemannianGeometry/#logdet-zero-1",
    "page": "Intro to Riemannian Geometry",
    "title": "logdet zero",
    "category": "section",
    "text": "The logdet zero divergence, also known as the square of the Bhattacharyya divergence (Mohaker, 2013)🎓, Stein divergence (Harandi et al., 2016)🎓, symmetrized Jensen divergence, the S-divergence (Sra, 2016)🎓 or the log determinant α-divergence (with α=0, Chebby and Moakher, 2012 🎓) is a Jensen-Bregman symmetric divergence enjoying all three invariances we have listed.Its square root has been shown to be a distance (Sra, 2016)🎓. It behaves very similarly to the Fisher metric at short distances (Moakher, 2012; Sra, 2016; Cichocki et al., 2015; Harandi et al., 2016) 🎓 and the mean of two matrices in P is the same as the Fisher mean  (Harandi et al., 2016) 🎓. Thus, it has often been used instead of the Fisher metric because it allows more efficient calculations. In fact, the calculation of this distance requires only three Cholesky decompositions, whereas the computation of the Fisher distance involves extracting generalized eigenvalues.distance² to I distance²\ntextrmlogdetfrac12(P+I)-frac12textrmlogdet(P) textrmlogdetfrac12(P+Q)-frac12textrmlogdet(PQ)geodesic: we use the Fréchet mean with appropriate weights.Fréchet mean: the solution is the unique positive definite matrix G satisfyingsum_i=1^kw_ibig(frac12P_i+frac12Gbig)^-1=G^-1.For estimating it PosDefManifold implements the fixed-point iterations (Moakher, 2012, p315)🎓:G  frack2big(sum_i=1^kw_i(P_i+G)^-1big)^-1.The logdet zero divergence between P and Q can also be written as the log-determinant of their arithmetic mean minus the log-determinant of their geometric mean (Moakher, 2012)🎓, which thus defines a possible extension to matrices of the useful concept of Wiener entropy."
},

{
    "location": "introToRiemannianGeometry/#logdet-α-1",
    "page": "Intro to Riemannian Geometry",
    "title": "logdet α",
    "category": "section",
    "text": "The log determinant α-divergence family for α-11 (Chebby and Moakher, 2012)🎓 allowsthe logdet zero mean for α=0,\nthe left Kullback-Leibler mean for α=-1 (which is the harmonic mean)\nthe right Kullback-Leibler mean for α=1 (which is the arithmetic mean).We do not consider the left and right Kullback-Leibler divergences because the related means are trivially the arithmetic and harmonic one (Moakher, 2012). As per the symmetrized Kullback-Leibler divergence, this is known as Jeffrey divergence and will be considered next. The log determinant α-divergence family of means is not implemented in PosDefManifold (besides the special cases α=(-1 0 1), since the family of power means are implemented."
},

{
    "location": "introToRiemannianGeometry/#Jeffrey-1",
    "page": "Intro to Riemannian Geometry",
    "title": "Jeffrey",
    "category": "section",
    "text": "This is a Jensen-Bregman symmetric divergence, also known as the symmetrized Kullback-Leibler divergence (see logdet α) (Faraki et al., 2015)🎓. It enjoyes all three invariances we have listed.distance² to I distance²\nfrac12textrmtr big(P+P^-1big)-n frac12textrmtr(Q^-1P+P^-1Q)-ngeodesic: we use the Fréchet mean with appropriate weights.Fréchet mean: A^12big(A^-12HA^-12big)^12A^12, where A is the arithmetic mean (see Euclidean metric) and H is the harmonic mean (see inverse Euclidean metric). Thus, the weighted Fréchet mean is the geometric mean (see Fisher metric) of the arithmetic and harmonic mean (Moakher, 2012)🎓.Note that this is the geometric mean only for k=2, that is, for scalars, but not in general for matrices, the geometric mean is the geometric mean of the arithmetic mean and harmonic mean (the only metric inducing the geometric mean in general is the Fisher mean)."
},

{
    "location": "introToRiemannianGeometry/#Von-Neumann-1",
    "page": "Intro to Riemannian Geometry",
    "title": "Von Neumann",
    "category": "section",
    "text": "The Von Neumann divergence is a Jensen-Bregman symmetric divergence (Sra, 2016; Taghia et al., 2019)🎓. It enjoyes only the rotation invariance.distance² to I distance²\nfrac12textrmtr(PtextrmlogP-textrmlogP) frac12textrmtrbig(P(textrmlogP-textrmlogQ)+Q(textrmlogQ-textrmlogP)big)The geodesic and weighted Fréchet mean for this metric are not available."
},

{
    "location": "introToRiemannianGeometry/#Wasserstein-1",
    "page": "Intro to Riemannian Geometry",
    "title": "Wasserstein",
    "category": "section",
    "text": "This is an extension to matrices of the Hellinger divergence for vectors and is also known as the Bures divergence in quantum physics, where it is applied on density matrices (unit trace positive-definite matrices). It enjoyes only the rotation invariance. Endowed with the Wasserstein metric the manifold P has a Riemannian geometry of nonnegative curvature. See ( Bhatia et al., 2019a; Bhatia et al., 2019b)🎓.distance² to I distance²\ntextrmtr(P+I)-2textrmtr(P^12) textrmtr(P+Q) -2textrmtrbig(P^12QP^12big)^12geodesic\nb^2P+a^2Q +abbig(PQ)^12 +(QP)^12bigThe quantity textrmtrbig(P^12QP^12big)^12 is known in quantum physics as the fidelity of P and  Q when those are density matrices (unit-trace positive definite matrices).Fréchet mean: the solution is the unique positive definite matrix G satisfying (Agueh and Carlier, 2011) 🎓G=sum_i=1^kw_ibig( G^12  P_i G^12big)^12.For estimating it, PosDefManifold implements the fixed-point algorithm of Álvarez-Esteban et al. (2016)🎓, giving iterations:G  G^-12 big(sum_i=1^k w_i(G^12P_i G^12)^12big)^2 G^-12In the special case when the matrices all pair-wise commute, the Wasserstein mean is equal to the instance of power means and generalized means with p=frac12 (Bhatia, Jain and Lim, 2019b)🎓, that is, to the modified Bhattacharyya mean.In the special case k=2 and equal weight the mean is W=frac14big(P+Q+(PQ) ^12+(QP)^12big)."
},

{
    "location": "introToRiemannianGeometry/#-1",
    "page": "Intro to Riemannian Geometry",
    "title": "🎓",
    "category": "section",
    "text": "ReferencesM. Agueh, G. Carlier (2011) Barycenters in the Wasserstein space, SIAM J. Mat. Anal. Appl. 43, 904-924.P. C. Álvarez-Esteban, E. del Barrio, J.A. Cuesta-Albertos, C. Matrána (2016) A fixed-point approach to barycenters in Wasserstein space, Journal of Mathematical Analysis and Applications, 441(2), 744-762.T. Ando, C.-K. Li, R. Mathias (2004) Geometric means, Linear Algebra and its Applications, 385(1), 305-334.V. Arsigny, P. Fillard, X. Pennec, N. Ayache (2007) Geometric means in a novel vector space structure on symmetric positive-definite matrices, SIAM journal on matrix analysis and applications, 29(1), 328-347.A. Barachant, S. Bonnet, M. Congedo, C. Jutten (2012) Multi-class Brain Computer Interface Classification by Riemannian Geometry, IEEE Transactions on Biomedical Engineering, 59(4), 920-928.A. Barachant, S. Bonnet, M. Congedo, C. Jutten (2013) Classification of covariance matrices using a Riemannian-based kernel for BCI applications, Neurocomputing, 112, 172-178.R. Bhatia (2007) Positive Definite Matrices. Princeton University press.R. Bhatia, M. Congedo (2019) Procrustes problems in manifolds of positive definite matrices Linear Algebra and its Applications, 563, 440-445.R. Bhatia, S. Gaubert, T. Jain (2019) Matrix versions of the Hellinger distance, arXiv:1901.01378.R. Bhatia, J. Holbrook (2006) Riemannian geometry and matrix geometric means, Linear Algebra and its Applications, 413 (2-3), 594-618.R. Bhatia, T. Jain (2010) Approximation problems in the Riemannian metric on positive definite matrices, Ann. Funct. Anal., 5(2), 118-126.R. Bhatia, T. Jain,Y. Lim (2019a) Inequalities for the Wasserstein mean of positive definite matrices, Linear Algebra and its Applications, in press.R. Bhatia, T. Jain, Y. Lim (2019b) On the Bures-Wasserstein distance between positive definite matrices Expositiones Mathematicae, in press.Z. Chebbi, M. Moakher (2012) Means of Hermitian positive-definite matrices based on the log-determinant α-divergence function, Linear Algebra and its Applications, 436(7), 1872-1889.A. Cichocki, S. Cruces, S-I- Amari (2015) Log-Determinant Divergences Revisited: Alpha-Beta and Gamma Log-Det Divergences, Entropy, 17(5), 2988-3034.M. Congedo, B. Afsari, A. Barachant, M Moakher (2015) Approximate Joint Diagonalization and Geometric Mean of Symmetric Positive Definite Matrices, PLoS ONE 10(4): e0121423.M. Congedo, A. Barachant, R. Bhatia R (2017a) Riemannian Geometry for EEG-based Brain-Computer Interfaces; a Primer and a Review, Brain-Computer Interfaces, 4(3), 155-174.M. Congedo, A. Barachant, E. Kharati Koopaei (2017b) Fixed Point Algorithms for Estimating Power Means of Positive Definite Matrices, IEEE Transactions on Signal Processing, 65(9), 2211-2220.X. Dai, S. Khamis, Y. Zhang, L.S. Davis (2016) Parameterizing region covariance: an efficient way to apply sparse codes on second order statistics, arXiv:1602.02822.M. Faraki, M. Harandi, F. Porikli (2015) More About VLAD: A Leap from Euclidean to Riemannian Manifolds, IEEE Conference on Computer Vision and Pattern Recognition (CVPR), Boston.W. Förstner, B. Moonen (1999) A metric for covariance matrices, In Krumm K and Schwarze VS eds. Qho vadis geodesia...?, number 1999.6 in tech. report of the Dep. Of Geodesy and Geoinformatics, p.113–128, Stuttgart University.M.T. Harandi, R. Hartley, B. Lovell, C. Sanderson (2016) Sparse coding on symmetric positive definite manifolds using bregman divergences, IEEE transactions on neural networks and learning systems, 27 (6), 1294-1306.Y. Lim, M. Pálfia (2012) Matrix power means and the Karcher mean,   Journal of Functional Analysis, 262(4), 1498-1514.Z. Lin (2019) Riemannian Geometry of Symmetric Positive Definite Matrices via Cholesky Decomposition. In press.M. Moakher (2005) A Differential Geometric Approach to the Geometric Mean of Symmetric Positive-Definite Matrices, SIAM Journal on Matrix Analysis and Applications, 26(3), 735-747.M. Moakher (2012) Divergence measures and means of symmetric positive-definite matrices, in D.H Lailaw and A. Vilanova (Eds) \"New Developments in the Visualization and Processing of Tensor Fields\", Springer, Berlin.X. Pennec, P. Fillard, N. Ayache (2006) A Riemannian Framework for Tensor Computing, International Journal of Computer Vision, 66(1), 41-66.P.L.C. Rodrigues, M. Congedo, C Jutten (2018) Multivariate Time-Series Analysis Via Manifold Learning, in Proc. of the the IEEE Statistical Signal Processing Workshop (SSP 2018), Fribourg-en-Brisgau, Germany.S. Sra (2016) Positive definite matrices and the S-divergence, Proc. Amer. Math. Soc., 144, 2787-2797.J. Taghia, M. Bånkestad, F. Lindsten, T.B. Schön (2019) Constructing the Matrix Multilayer Perceptron and its Application to the VAE, arXiv:1902.01182v1S. Umeyama (1988) An Eigendecomposition Approach to Weighted Graph Matching Problems, IEEE Trans. Pattern. Anal. Mach. Intell., 10(5), 695-703."
},

{
    "location": "MainModule/#",
    "page": "MainModule (PosDefManifold.jl)",
    "title": "MainModule (PosDefManifold.jl)",
    "category": "page",
    "text": ""
},

{
    "location": "MainModule/#MainModule-(PosDefManifold.jl)-1",
    "page": "MainModule (PosDefManifold.jl)",
    "title": "MainModule (PosDefManifold.jl)",
    "category": "section",
    "text": "This is the main unit containing the PosDefManifold module.It uses the following standard Julia packages:using\nLinear Algebra\nStatisticsExamples in some units of PosDefManifold also uses the Plots package.The main module does not contains functions, but it declares all constant, types and aliases of Julia functions used in all units.Contents\nconstants\ntypes\naliases\ntips & tricks"
},

{
    "location": "MainModule/#constants-1",
    "page": "MainModule (PosDefManifold.jl)",
    "title": "constants",
    "category": "section",
    "text": "constant value numeric value\nsqrt2 √2 1.4142135623730951\ninvsqrt2 1/√2 0.7071067811865475\nminpos 1e-15 0.000000000000001\nmaxpos 1e15 100000000000000"
},

{
    "location": "MainModule/#types-1",
    "page": "MainModule (PosDefManifold.jl)",
    "title": "types",
    "category": "section",
    "text": ""
},

{
    "location": "MainModule/#Metric::Enumerated-type-1",
    "page": "MainModule (PosDefManifold.jl)",
    "title": "Metric::Enumerated type",
    "category": "section",
    "text": "@enum Metric begin\n  Euclidean    =1\n  invEuclidean =2\n  ChoEuclidean =3\n  logEuclidean =4\n  LogCholesky  =5\n  Fisher       =6 # default metric\n  logdet0      =7\n  Jeffrey      =8\n  VonNeumann   =9\n  Wasserstein  =10\nendRiemannian manipulations are defined for a given metric (see metrics).  An instance for this type is requested as an argument in many functions  contained in the riemannianGeometry.jl unit in order to specify  the metric, unless the default metric (Fisher) is sought. ## Example\n # generate a 15x15 symmetric positive definite matrix\n P=randP(15)              \n # distance from P to the identity matrix according to the logdet0 metric\n d=distance(P, logdet0)  If you want to work consistently with a specific metric,  you may want to declare in your script a global variable such asglobal metric=logdet0  or  global metric=Metric(Int(logdet0)),and then pass metric as argument in all your computations,  e.g., referring to the above example,d=distance(P, metric).To know what is the current metric, get it as a string as:s=string(metric)"
},

{
    "location": "MainModule/#RealOrComplex-type-1",
    "page": "MainModule (PosDefManifold.jl)",
    "title": "RealOrComplex type",
    "category": "section",
    "text": "RealOrComplex=Union{Real, Complex} is the Union of Real and Complex Types."
},

{
    "location": "MainModule/#aliases-1",
    "page": "MainModule (PosDefManifold.jl)",
    "title": "aliases",
    "category": "section",
    "text": "alias Julia function in Package tab-completition REPL support\n𝚺 sum Base \\bfSigma ⛔\n𝛍 mean Statistics \\bfmu ⛔\nℂ ComplexF64 Base \\bbC ✓\n⋱ Diagonal LinearAlgebra \\ddots ✓\nℍ Hermitian LinearAlgebra \\bbH ✓All packages above are built-in julia packages."
},

{
    "location": "MainModule/#tips-and-tricks-1",
    "page": "MainModule (PosDefManifold.jl)",
    "title": "tips & tricks",
    "category": "section",
    "text": ""
},

{
    "location": "MainModule/#typecasting-matrices-1",
    "page": "MainModule (PosDefManifold.jl)",
    "title": "typecasting matrices",
    "category": "section",
    "text": "Several functions in PosDefManifold implement multiple dispatch and can handle    several kinds of matrices as input, however the core functions for manipulating    objects on the Riemannian manifold of positive definite matrices act by definition  on positive definite matrices only.  Those matrices must therefore be either  symmetric positive definite (real) or Hermitian (complex).  Such matrices are identified in as being of the Hermitiantype, using the standard LinearAlgebra package.  The alias ℍ is used consistently in the code (see aliases).  If the input is not flagged, the functions restricting the input to  positive definite matrices will give an error.Examplejulia> using LinearAlgebra\n\njulia> f(S::Hermitian)=S*S\'\nf (generic function with 1 method)\n\njulia> A=randn(3, 3)\n3×3 Array{Float64,2}:\n -0.67407  -0.344258    0.203714\n -1.06551  -0.0233796   0.975465\n -1.04727  -1.19807    -0.0219121\n\njulia> H=A*A\' # although SPD, H is not automatically flagged as Hermitian\n3×3 Array{Float64,2}:\n 0.614384  0.924991  1.11391\n 0.924991  2.08738   1.12251\n 1.11391   1.12251   2.53263\n\njulia> f(H)\nERROR: MethodError: no method matching f(::Array{Float64,2})\nClosest candidates are:\n  f(::Hermitian) at none:1If you construct a positive definite matrix and it is not flagged,  you can do so simply by typecasting it, that is, passing as argument to the  functions Hermitian(P) instead of just P. The ℍ alias can be  used for short, i.e., ℍ(P). Continuing the example above:julia> f(ℍ(H))  # this way it works, equivalent to f(Hermitian(H))\n3×3 Array{Float64,2}:\n 2.47388  3.74948  4.54381\n 3.74948  6.4728   6.21635\n 4.54381  6.21635  8.91504Finally, other functions act on generic matrices (of type Matrix).  To those functions you can pass any matrix.  However, keep in mind that the functions writing on the argument matrix such as  normalizeCol! will give an error if you pass an Hermitian matrix,  since Julia does not allow writing on non-diagonal elements of those matrices.  In this case typecast it in another object using the Matrix type;  Suppose H is Hermitian, you would use for example:julia> X=Matrix(H)\njulia> normalizeCol!(X, 1)\njulia> norm(X[:, 1])\n1.0Another example when typecasting is useful: the gram function  takes a Matrix type as argument (since X is expected to be a data matrix),  like inH=gram(X)The following will not work though: H=gram(X\')since X\' is an Adjoint type. The problem is fixed by typecasting the  adjoint matrix, such asH=gram(Matrix(X\'))Another example: here is how to get an Hermitian matrix out of the  diagonal part of an Hermitian matrix H:Hermitian(Matrix(Diagonal(H)))"
},

{
    "location": "riemannianGeometry/#",
    "page": "riemannianGeometry.jl",
    "title": "riemannianGeometry.jl",
    "category": "page",
    "text": ""
},

{
    "location": "riemannianGeometry/#riemannianGeometry.jl-1",
    "page": "riemannianGeometry.jl",
    "title": "riemannianGeometry.jl",
    "category": "section",
    "text": "This is the fundamental unit of PosDefManifold. It contains functions for manipulating points in the Riemannian manifold of Symmetric Positive Definite (SPD) or Hermitian Positive Definite (HPD) matrices. In Julia those are Hermitian matrices, see typecasting matrices.The functions are divided in six categories:Category Output\n1. Geodesic equations interpolation, extrapolation,...\n2. Distances length of geodesics\n3. Graphs and Laplacians for spectral embedding, eigenmaps, system dynamics,...\n4. Means mid-points of geodesics, centers of mass of several points\n5. Tangent Space operations maps from the manifold to the tangent space and viceversa\n6. Procrustes problems for data matching, transfer learning,...⋅"
},

{
    "location": "riemannianGeometry/#Geodesic-equations-1",
    "page": "riemannianGeometry.jl",
    "title": "Geodesic equations",
    "category": "section",
    "text": "geodesic"
},

{
    "location": "riemannianGeometry/#Distances-1",
    "page": "riemannianGeometry.jl",
    "title": "Distances",
    "category": "section",
    "text": "distanceSqr\ndistance"
},

{
    "location": "riemannianGeometry/#Graphs-and-Laplacians-1",
    "page": "riemannianGeometry.jl",
    "title": "Graphs and Laplacians",
    "category": "section",
    "text": "distanceSqrMat\ndistanceMatrix\nlaplacian\nlaplacianEigenMaps\nspectralEmbedding"
},

{
    "location": "riemannianGeometry/#Means-1",
    "page": "riemannianGeometry.jl",
    "title": "Means",
    "category": "section",
    "text": "generalizedMean\nlogdet0Mean\nwasMean\npowerMean\nmeanP"
},

{
    "location": "riemannianGeometry/#Tangent-Space-operations-1",
    "page": "riemannianGeometry.jl",
    "title": "Tangent Space operations",
    "category": "section",
    "text": "logMap\nexpMap\nvecP\nmatP"
},

{
    "location": "riemannianGeometry/#Procrustes-problems-1",
    "page": "riemannianGeometry.jl",
    "title": "Procrustes problems",
    "category": "section",
    "text": "procrustes"
},

{
    "location": "linearAlgebra/#",
    "page": "linearAlgebra.jl",
    "title": "linearAlgebra.jl",
    "category": "page",
    "text": ""
},

{
    "location": "linearAlgebra/#linearAlgebra.jl-1",
    "page": "linearAlgebra.jl",
    "title": "linearAlgebra.jl",
    "category": "section",
    "text": "This unit contains linear algebra functions useful in relation to the Riemannian  geometry of the manifold of Symmetric Positive Definite (SPD) or  Hermitian Positive Definite (HPD) matrices. In Julia those are Hermitian matrices, see typecasting matrices.In general they take a matrix as input (some may take other arrays as input) and are divided in seven categories depending on what kind of functions thay are and what they give as output:Category Output\n1. Matrix normalizations matrix\n2. Boolean functions of matrices matrix\n3. Scalar functions of matrices scalar\n4. Diagonal functions of matrices diagonal matrix\n5. Unitary functions of matrices orthogonal/unitary matrix\n6. Matrix function of matrices matrix\n7. Spectral decompositions of positive matrices spectral function of input\n8. Decompositions involving triangular matrices triangular matrix⋅"
},

{
    "location": "linearAlgebra/#Matrix-normalizations-1",
    "page": "linearAlgebra.jl",
    "title": "Matrix normalizations",
    "category": "section",
    "text": "Function Description\ndet1 Normalize the determinant\ntr1 Normalize the trace\nnormalizeCol! Normalize one or more columns⋅det1\ntr1\nnormalizeCol!"
},

{
    "location": "linearAlgebra/#Boolean-functions-of-matrices-1",
    "page": "linearAlgebra.jl",
    "title": "Boolean functions of matrices",
    "category": "section",
    "text": "Function Description\nispos Check whether the argument is comprised of all positive elementsispos"
},

{
    "location": "linearAlgebra/#Scalar-functions-of-matrices-1",
    "page": "linearAlgebra.jl",
    "title": "Scalar functions of matrices",
    "category": "section",
    "text": "Function Description\ncolProd Sum of products of the elements in two columns\nsumOfSqr Sum of squares of all elements or of specified columns\nsumOfSqrDiag Sum of squares of the diagonal elements\ncolNorm Eucliden norm of a column\nsumOfSqrTril Sum of squares of the lower triangle elements up to a given underdiagonal\nfidelity (Quantum) Fidelity of two positive matrices⋅colProd\nsumOfSqr\nsumOfSqrDiag\ncolNorm\nsumOfSqrTril\nfidelity"
},

{
    "location": "linearAlgebra/#Diagonal-functions-of-matrices-1",
    "page": "linearAlgebra.jl",
    "title": "Diagonal functions of matrices",
    "category": "section",
    "text": "Function Description\nfDiagonal Elemen-wise functions of matrix diagonals⋅fDiagonal"
},

{
    "location": "linearAlgebra/#Unitary-functions-of-matrices-1",
    "page": "linearAlgebra.jl",
    "title": "Unitary functions of matrices",
    "category": "section",
    "text": "Function Description\nmgs Modified Gram-Schmidt orthogonalization⋅mgs"
},

{
    "location": "linearAlgebra/#Matrix-function-of-matrices-1",
    "page": "linearAlgebra.jl",
    "title": "Matrix function of matrices",
    "category": "section",
    "text": "Function Description\nnone for now ipse lorem...⋅"
},

{
    "location": "linearAlgebra/#Spectral-decompositions-of-positive-matrices-1",
    "page": "linearAlgebra.jl",
    "title": "Spectral decompositions of positive matrices",
    "category": "section",
    "text": "Function Description\nevd Eigenvalue-Eigenvector decomposition of a matrix in UΛU=P form\nspectralFunctions Mother function for creating spectral functions of eigenvalues\npow Power of a positive matrix for any number of exponents in one pass\ninvsqrt Principal square root inverse (whitening) of a positive matrix\nsqr Square of a positive matrix\npowerIterations, powIter Power method for estimating any number of eigenvectors and associated eigenvalues\nchoL Lower triangula factor of Cholesky decomposition⋅evd\nspectralFunctions\npow\ninvsqrt\nsqr\npowerIterations"
},

{
    "location": "linearAlgebra/#Decompositions-involving-triangular-matrices-1",
    "page": "linearAlgebra.jl",
    "title": "Decompositions involving triangular matrices",
    "category": "section",
    "text": "choL"
},

{
    "location": "signalProcessing/#",
    "page": "signalProcessing.jl",
    "title": "signalProcessing.jl",
    "category": "page",
    "text": ""
},

{
    "location": "signalProcessing/#signalProcessing.jl-1",
    "page": "signalProcessing.jl",
    "title": "signalProcessing.jl",
    "category": "section",
    "text": "This unit contains miscellaneous signal processing functions useful in relation to the Riemannian geometry of the manifold of Symmetric Positive Definite (SPD) or Hermitian Positive Definite (HPD) matrices. In Julia those are Hermitian matrices, see typecasting matrices.Function Description\nrandChi², randχ² Generate a random variable distributed as a chi-squared\nrandEigvals, randλ Generate a random vectors of real positive eigenvalues\nrandEigvalsMat, randΛ Generate a random diagonal matrix of real positive eigenvalues\nrandUnitaryMat, randU Generate a random orthogonal or unitary matrix\nrandPosDefMat, randP Generate one or an array of random positive definite matrices\nregularize! Regularize an array of positive definite matrices\ngram Gram matrix of a matrix\ntrade trace and determinant of a matrix as a 2-tuple⋅randChi²\nrandEigvals\nrandEigvalsMat\nrandUnitaryMat\nrandPosDefMat\nregularize!\ngram\ntrade"
},

{
    "location": "test/#",
    "page": "test.jl",
    "title": "test.jl",
    "category": "page",
    "text": ""
},

{
    "location": "test/#test.jl-1",
    "page": "test.jl",
    "title": "test.jl",
    "category": "section",
    "text": "Most functions in PosDefManifold are tested, both for real and complex data input. This unit declares the function testall() that performs all tests.The first time you execute the test it will take some time. If a test fails you will be pointed to the function that failed."
},

]}

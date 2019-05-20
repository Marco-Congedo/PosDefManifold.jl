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
    "location": "#Requirements-1",
    "page": "PosDefManifold Documentation",
    "title": "Requirements",
    "category": "section",
    "text": "Julia version ≥ 1.0.3"
},

{
    "location": "#Installation-1",
    "page": "PosDefManifold Documentation",
    "title": "Installation",
    "category": "section",
    "text": "The package is not registered yet. Execute the following command in Julia\'s REPL:]add https://github.com/Marco-Congedo/PosDefManifold.jl"
},

{
    "location": "#Overview-1",
    "page": "PosDefManifold Documentation",
    "title": "Overview",
    "category": "section",
    "text": "(Image: Figure 1)Riemannian geometry studies smooth manifolds, multi-dimensional curved spaces with peculiar geometries endowed with non-Euclidean metrics. In these spaces Riemannian geometry allows the definition of angles, geodesics (shortest path between two points), distances between points, centers of mass of several points, etc.In this package we are concerned with the manifold P of positive definite matrices, either symmetric positive definite or Hermitian positive definite.In several fields of research such as computer vision and brain-computer interface, treating data in the P manifold has allowed the introduction of machine learning approaches with remarkable characteristics, such as simplicity of use, excellent classification accuracy, as demonstrated by the winning score obtained in six international data classification competitions, and the ability to operate transfer learning (Congedo et al., 2017)🎓).For a formal introduction to the P manifold the reader is referred to the monography written by Bhatia (2007)🎓.For an introduction to Riemannian geometry and an overview of mathematical tools implemented in this package, see Intro to Riemannian Geometry in this documentation.For starting using this package, browse the code units listed here below and execute the many code examples you will find therein. The core functions are contained in unit riemannianGeometry.jl."
},

{
    "location": "#Code-units-1",
    "page": "PosDefManifold Documentation",
    "title": "Code units",
    "category": "section",
    "text": "PosDefManifold includes six code units (.jl files):Unit Description\nMainModule (PosDefManifold.jl) Main module, constants, types, aliases, tips & tricks\nriemannianGeometry.jl The fundamental unit collecting all functions acting on the P manifold\nlinearAlgebra.jl Collection of linear algebra routines\nsignalProcessing.jl Collection of signal processing routines\nclassification.jl Collection of classification routines\ntest.jl Unit performing all tests"
},

{
    "location": "#Contents-1",
    "page": "PosDefManifold Documentation",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [       \"index.md\",\n                \"introToRiemannianGeometry.md\",\n                \"MainModule.md\",\n                \"riemannianGeometry.md\",\n                \"linearAlgebra.md\",\n                \"signalProcessing.md\",\n                \"classification.md\",\n                \"test.md\"]\nDepth = 1"
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
    "text": "The study of appropriate distance measures for positive definite matrices has recently grown very fast, driven by practical problems in radar data processing, image processing, computer vision, shape analysis, medical imaging (especially diffusion MRI and Brain-Computer Interface), sensor networks, elasticity, mechanics, numerical analysis and machine learning (e.g., see references in Congedo et al., 2017a)🎓.In many applications the observed data can be conveniently summarized by positive definite matrices, which are either symmetric positive definite (SPD: real) or Hermitian Positive Definite (HPD: complex). For example, those may be some form of the data covariance matrix in the time, frequency or time-frequency domain, or autocorrelation matrices, kernels, slices of tensors, density matrices, elements of a search space, etc. Positive definite matrices are naturally treated as points on a smooth Riemannian manifold allowing useful operations such as interpolation, smoothing, filtering, approximation, averaging, signal detection and classification. Such operations are the object of the present PosDefManifold library.More formally, this Julia library treats operations on the metric space (P δ^2) of n・n positive definite matrices endowed with a distance or symmetric divergence δ(P x P)0 . Several matrix distances or matrix divergences δ are considered. Using some of them, the most important one being the Fisher metric, we define a Riemannian manifold. In mathematics, this is the subject of Riemannian geometry and information geometry.Note that throughout this library the word \'metric\' is used loosely for referring to the actual Riemannian metric on the tangent space and to the resulting distance or to general symmetric divergence acting on P, regardless the fact that we are dealing with a metric in the strict sense and that it induces or not a Riemannian geometry in P. This is done for convenience of exposition, since in practice those \'metrics\' in PosDefManifold may be used interchangeably."
},

{
    "location": "introToRiemannianGeometry/#Riemannian-manifolds-1",
    "page": "Intro to Riemannian Geometry",
    "title": "Riemannian manifolds",
    "category": "section",
    "text": "Here are some important definitions:A smooth manifold in differential geometry is a topological space that is locally similar to the Euclidean space and has a globally defined differential structure.The tangent space at point G is the vector space containing the tangent vectors to all curves on the manifold passing through G (Fig. 1).A smooth Riemannian manifold is equipped with an inner product on the tangent space (a Riemannian metric) defined at each point and varying smoothly from point to point. For manifold P the tangent space is the space of symmetric or Hermitian matrices.Thus, a Riemannian metric turns the metric space (P δ^2) into a Riemannian manifold. This is the case, for example, of the Fisher metric, which has a fundamental role in the manifolds of positive definite matrices and of the Wasserstein metric, fundamental in optimal transport theory.(Image: Figure 1) Figure 1. Schematic illustration of the Riemannian manifold of positive definite matrices. Left: geodesic relying points P and Q passing through its-mid-point (mean) G (green curve), tangent space at point G with tangent vectors to geodesic from G to P and from G to Q (blue arrowed lines) and distance δ(G Q). Right: the center of mass (also named mean) G of points P1P4 defined as the point minimizing the sum of the four squared distances δ2(G P_i), for i=14."
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
    "text": "The length of the geodesic (at constant velocity) between two points gives the distance δ(P Q).  The distance is always real, non-negative and equal to zero if and only if P=Q."
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
    "text": "An important characteristic of metrics is that they may induce invariance properties on the distance, which are in turn inherited by the mean.Let us denote shortly by P_i the set P_1P_k, where i=1k  and by GP_i the Fréchet mean of the set (in this section we drop the weights here for keeping the notation short). The most important invariance properties are:invariance effect on distance δ(PQ) effect on mean GP_i\nrotation δ(PQ)=δ(U^*PUU^*QU) GU^*P_iU=U^*GP_iU\naffinity δ(PQ)=δ(B^*PBB^*QB) GB^*P_iB=B^*GP_iB\ninversion δ(PQ)=δ(P^-1Q^-1) GP_i^-1=G^-1P_ifor any unitary U unitary and non-singular B.The affine invariance implies the rotation invariance and is also named congruence invariance."
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
    "text": "The Fisher metric, also known as affine-invariant, natural and Fisher-Rao metric, among others names, has a paramount importance for the P manifold, standing out as the natural choice both from the perspective of differential geometry and information geometry. Endowed with the Fisher metric the manifold P is Riemannian, has nonpositive curvature and is symmetric. This metric verifies all three invariances we have considered.Distance² to I Distance²\ntextrmlog(P)^2 textrmlog(P^-12QP^-12)^2geodesic\nP^12 big(P^-12 Q P^-12big)^a P^12Fréchet mean: it does not have a closed-form solution in general. The solution is the unique positive definite matrix G satisfying (Bhatia and Holbrook, 2006; Moakher, 2005).🎓sum_i=1^kw_itextrmlogbig(G^-12 P_i G^-12big)=0For estimating it, PosDefManifold implements the well-known gradient descent algorithm, resulting in iterations:G G^12textrmexpbig(sum_i=1^kw_itextrmlog(G^-12 P_i G^-12)big)G^12Alternatively, and more efficiently, one can ask for an approximate solution invoking the MPM algorithm (Congedo et al., 2017b)🎓, which is also implemented (in order to estimate the geometric mean use function powerMean with parameter p=0 or with a very small value of p).This mean is known under many different names (Fisher, Rao, Fisher-Rao, Pusz-Woronowicz, Cartan, Fréchet, Karcher, geometric....).   The ‘centrality’ of this mean among a wide family of divergence-based means can be appreciated in Fig. 4 of Cichocki et al. (2015)🎓.The geometric mean G of two matrices P and Q is denoted shortly as PtextrmQ. Currently it is an object of intense study because of its interesting mathematical properties. For instance,it is the unique solution to Riccati equation GQ^-1G=P\nit is equal to F^-*D_1^12D_2^12F^-1 for whatever joint diagonalizer F of P and Q, i.e., for whatever matrix F satisfying F^*PF=D_1 and F^*QF=D_2, with D_1, D_1 non-singular diagonal matrices (Congedo et al., 2015)🎓.\nit enjoys all 10 properties of means postulated in the seminal work of Ando et al. (2010)🎓.When P and Q commutes, the Fisher mean of two matrices reduces to P^12Q^12, which indeed in this case is the log-Euclidean mean frac12textrmlogP + frac12textrmlogQ.The Fisher geodesic equation is usually denoted Ptextrm_aQ. Note that Itextrm_aP=P^a and Ptextrm_aI=P^b, where b=1-a.Fisher geodesic equation verifies Ptextrm_aQ=Qtextrm_bP and (Ptextrm_aQ)^-1=P^-1textrm_aQ^-1.An interesting property of the Fisher metric is that using its geodesic equation we can extrapolate positive matrices, always remaining in P. That is, using any real value of a :with 0  a  1 we move toward Q		(attraction),\nwith a  1 we move over and beyond Q	(extrapolation) and\nwith a 0 we move back away from Q 	(repulsion).Something similar can be done using the log Cholesky metric as well."
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
    "text": "If matrices P_1 P_2P_k all pair-wise commute, the special case p=frac12 yields the following instance of power means (and of generalized means):big(sum_i=1^kw_iP_i^12big)^12.This mean has been proposed  in a different context by Moakher (2012)🎓 as a modified Bhattacharyya mean, since it is a modification of the Bhattacharyya mean we will encounter next under the name logdet zero. It is worth noting that in commuting algebra Moakher’s mean also corresponds to the mean obtained with the Wasserstein metric."
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
    "text": "This is the main unit containing the PosDefManifold module.It uses the following standard Julia packages:using\nLinearAlgebra\nStatisticsExamples in some units of PosDefManifold also uses the Plots package. Take a look at this tutorial for an introduction to data plotting with Julia.The main module does not contains functions, but it declares all constant, types and aliases of Julia functions used in all units.Contents\nconstants\naliases\ntypes\ntips & tricks"
},

{
    "location": "MainModule/#constants-1",
    "page": "MainModule (PosDefManifold.jl)",
    "title": "constants",
    "category": "section",
    "text": "constant value numeric value\nsqrt2 √2 1.4142135623730951\ninvsqrt2 1/√2 0.7071067811865475\nmaxpos 1e15 100000000000000"
},

{
    "location": "MainModule/#aliases-1",
    "page": "MainModule (PosDefManifold.jl)",
    "title": "aliases",
    "category": "section",
    "text": "alias Julia function in Package tab-completition REPL support\n𝚺 sum Base \\bfSigma ⛔\n𝛍 mean Statistics \\bfmu ⛔\n𝕄 Matrix Base \\bbM ⛔\n𝔻 Diagonal LinearAlgebra \\bbD ⛔\nℍ Hermitian LinearAlgebra \\bbH ✓\n𝕃 LowerTriangular LinearAlgebra \\bbH ⛔All packages above are built-in julia packages."
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
    "text": "@enum Metric begin\n  Euclidean    =1\n  invEuclidean =2\n  ChoEuclidean =3\n  logEuclidean =4\n  LogCholesky  =5\n  Fisher       =6\n  logdet0      =7\n  Jeffrey      =8\n  VonNeumann   =9\n  Wasserstein  =10\nendRiemannian manipulations are defined for a given metric (see metrics).  An instance for this type is requested as an argument in many functions  contained in the riemannianGeometry.jl unit in order to specify  the metric. ## Example\n # generate a 15x15 symmetric positive definite matrix\n P=randP(15)              \n # distance from P to the identity matrix according to the logdet0 metric\n d=distance(logdet0, P)  If you want to work consistently with a specific metric,  you may want to declare in your script a global variable such asglobal metric=logdet0  or  global metric=Metric(Int(logdet0)),and then pass metric as argument in all your computations,  e.g., referring to the above example,d=distance(metric, P).To know what is the current metric, you can get it as a string using:s=string(metric)To see the list of metrics in type Metric use:instances(Metric)"
},

{
    "location": "MainModule/#Array-of-Matrices-types-1",
    "page": "MainModule (PosDefManifold.jl)",
    "title": "Array of Matrices types",
    "category": "section",
    "text": ""
},

{
    "location": "MainModule/#𝕄Vector-type-1",
    "page": "MainModule (PosDefManifold.jl)",
    "title": "𝕄Vector type",
    "category": "section",
    "text": "𝕄Vector=Vector{𝕄}This is a vector of general Matrix matrices, alias of MatrixVector.  Julia sees is at: Array{Array{T,2} where T,1}.See aliases for the 𝕄 symbol and typecasting matrices for the use of matrices in PosDefManifold.warning: Nota bene\nThis object is meant to hold matrices living in the same manifold, therefore it is assumed by all methods that all matrices it holds are of the same dimension.See dim, typeofMatrix𝕄Vector₂ type𝕄Vector₂=Vector{𝕄Vector}This is a vector of 𝕄Vector type objects,  i.e., a vector of vectors of matrices.  It is the alias of MatrixVector₂.  Julia sees it as: Array{Array{Array{T,2} where T,1},1}.warning: Nota bene\nThis object is meant to hold matrices living in the same manifold, therefore it is assumed by all methods that all matrices it holds are of the same dimension. However the several 𝕄Vector objects it holds do not need to have the same length.    See dim, typeofMatrix"
},

{
    "location": "MainModule/#𝔻Vector-type-1",
    "page": "MainModule (PosDefManifold.jl)",
    "title": "𝔻Vector type",
    "category": "section",
    "text": "𝔻Vector=Vector{𝔻}This is a vector of Diagonal matrices, alias of DiagonalVector.  Julia sees is at: Array{Diagonal,1}.See aliases for the 𝔻 symbol and typecasting matrices for the use of Diagonal matrices in PosDefManifold.warning: Nota bene\nThis object is meant to hold matrices living in the same manifold, therefore it is assumed by all methods that all matrices it holds  are of the same dimension.See dim, typeofMatrix𝔻Vector₂ type𝔻Vector₂=Vector{𝔻Vector}This is a vector of 𝔻Vector type objects,  i.e., a vector of vectors of Diagonal matrices.  It is the alias of DiagonalVector₂.  Julia sees it as: Array{Array{Diagonal,1},1}.warning: Nota bene\nThis object is meant to hold matrices living in the same manifold, therefore it is assumed by all methods that all matrices it holds are of the same dimension. However the several 𝔻Vector objects it holds do not need to have the same length.    See dim, typeofMatrix"
},

{
    "location": "MainModule/#𝕃Vector-type-1",
    "page": "MainModule (PosDefManifold.jl)",
    "title": "𝕃Vector type",
    "category": "section",
    "text": "𝕃Vector=Vector{𝕃}This is a vector of LowerTriangular matrices, alias of LowerTriangularVector.  Julia sees is at: Array{LowerTriangular,1}.See aliases for the 𝕃 symbol and typecasting matrices for the use of LowerTriangular matrices in PosDefManifold.warning: Nota bene\nThis object is meant to hold matrices living in the same manifold, therefore it is assumed by all methods that all matrices it holds are of the same dimension.See dim, typeofMatrix𝕃Vector₂ type𝕃Vector₂=Vector{𝕃Vector}This is a vector of 𝕃Vector type objects, i.e.,  a vector of vectors of LowerTriangular matrices.  It is the alias of LowerTriangularVector₂.  Julia sees it as: Array{Array{LowerTriangular,1},1}.warning: Nota bene\nThis object is meant to hold matrices living in the same manifold, therefore it is assumed by all methods that all matrices it holds are of the same dimension. However the several 𝕃Vector objects it holds do not need to have the same length.    See dim, typeofMatrix"
},

{
    "location": "MainModule/#ℍVector-type-1",
    "page": "MainModule (PosDefManifold.jl)",
    "title": "ℍVector type",
    "category": "section",
    "text": "ℍVector=Vector{ℍ}This is a vector of Hermitian matrices, alias of HermitianVector.  Julia sees is at: Array{Hermitian,1}.See aliases for the ℍ symbol and typecasting matrices for the use of Hermitian matrices in PosDefManifold.warning: Nota bene\nThis object is meant to hold matrices living in the same manifold, therefore it is assumed by all methods that all matrices it holds are of the same dimension.See dim, typeofMatrixℍVector₂ type`ℍVector₂=Vector{ℍVector}`This is a vector of ℍVector type  objects, i.e., a vector of vectors of Hermitian matrices.  It is the alias of HermitianVector₂.  Julia sees it as: Array{Array{Hermitian,1},1}.warning: Nota bene\nThis object is meant to hold matrices living in the same manifold, therefore it is assumed by all methods that all matrices it holds are of the same dimension. However the several ℍVector objects it holds do not need to have the same length.See dim, typeofMatrix"
},

{
    "location": "MainModule/#RealOrComplex-type-1",
    "page": "MainModule (PosDefManifold.jl)",
    "title": "RealOrComplex type",
    "category": "section",
    "text": "RealOrComplex=Union{Real, Complex}This is the Union of Real and Complex types."
},

{
    "location": "MainModule/#AnyMatrix-type-1",
    "page": "MainModule (PosDefManifold.jl)",
    "title": "AnyMatrix type",
    "category": "section",
    "text": "AnyMatrix=Union{𝔻{T}, 𝕃{T}, ℍ{T}, 𝕄{T}} where T<:RealOrComplexThis is the Union  of real or complex Diagonal, LowerTriangular, Hermitian and Matrix types. It is often used in the definition of functions.See aliases"
},

{
    "location": "MainModule/#AnyMatrixVector-type-1",
    "page": "MainModule (PosDefManifold.jl)",
    "title": "AnyMatrixVector type",
    "category": "section",
    "text": "AnyMatrixVector=Union{𝕄Vector, 𝔻Vector, 𝕃Vector, ℍVector}This is the Union of 𝕄Vector, 𝔻Vector, 𝕃Vector and ℍVector. It is often used in the definition of functions.  See Array of Matrices types.AnyMatrixVector₂ typeAnyMatrixVector₂=Union{𝕄Vector₂, 𝔻Vector₂, 𝕃Vector₂, ℍVector₂}This is the Union of 𝕄Vector₂, 𝔻Vector₂, 𝕃Vector₂, ℍVector₂. It is often used in the definition of functions. See Array of Matrices types."
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
    "text": "Several functions in PosDefManifold implement multiple dispatch and can handle several kinds of matrices as input, however the core functions for manipulating objects on the Riemannian manifold of positive definite matrices act by definition on positive definite matrices only.  Those matrices must therefore be either  symmetric positive definite (SPD, real) or Hermitian positive definite (HPD, complex).  Such matrices are uniformly identified in PosDefManifold as being of the Hermitian type, using the standard LinearAlgebra package.  The alias ℍ is used consistently in the code (see aliases).  If the input is not flagged as Hermitian, the functions restricting the input to positive definite matrices will not be accessible.Examplejulia> using LinearAlgebra\n\njulia> f(S::Hermitian)=S*S\'\nf (generic function with 1 method)\n\njulia> A=randn(3, 3)\n3×3 Array{Float64,2}:\n -0.67407  -0.344258    0.203714\n -1.06551  -0.0233796   0.975465\n -1.04727  -1.19807    -0.0219121\n\njulia> H=A*A\' # although SPD, H is not automatically flagged as Hermitian\n3×3 Array{Float64,2}:\n 0.614384  0.924991  1.11391\n 0.924991  2.08738   1.12251\n 1.11391   1.12251   2.53263\n\njulia> f(H)\nERROR: MethodError: no method matching f(::Array{Float64,2})\nClosest candidates are:\n  f(::Hermitian) at none:1If you construct a positive definite matrix and it is not flagged,  you can do so simply by typecasting it, that is, passing as argument to the  functions Hermitian(P) instead of just P. The ℍ alias can be  used for short, i.e., ℍ(P). Continuing the example above:julia> f(ℍ(H))  # this way it works, equivalent to f(Hermitian(H))\n3×3 Array{Float64,2}:\n 2.47388  3.74948  4.54381\n 3.74948  6.4728   6.21635\n 4.54381  6.21635  8.91504Be careful: Hermitian(P) will construct and Hermitian matrix from the argument.  If the matrix argument is not symmetric (if real) or Hermitian (if complex)  it will be made so by copying the transpose (if real) or complex conjugate  and transpose (if complex) of a triangular part into the other.  See Hermitian.If you want to construct an ℍVector type from, say, two Hermitian matrices P and Q, don\'t write A=[P, Q], but rather A=ℍVector([P, Q]). In fact,  the first is seen by Julia as2-element Array{Hermitian{Float64,Array{Float64,2}},1},while the latter as2-element Array{Hermitian,1},which is the type expected in all functions taking an ℍVector type as argument.Other functions act on generic matrices (of type Matrix). This is seen by Julia as Array{T,2} where T.  Keep in mind that the functions writing on the argument matrix such as  normalizeCol! will give an error if you pass an Hermitian matrix,  since Julia does not allow writing on non-diagonal elements of those matrices.  In this case typecast it in another object using the Matrix type;  suppose H is Hermitian, you would use for example:julia> X=Matrix(H)\njulia> normalizeCol!(X, 1)\njulia> norm(X[:, 1])\n1.0Some more examples:Typecasting Adjoint matrices:\nMatrix(X\')\nhere is how to get an Hermitian matrix out of thediagonal part of an Hermitian matrix H:Hermitian(Matrix(Diagonal(H)))here is how to get a LowerTriangular matrix out of anHermitian matrix H:LowerTriangular(Matrix(H))For example, you can use this to pass a full inter-distance matrix to the laplacian function to obtain the Laplacian matrix.A useful function is typeofMatrix. For example, the following line  typecasts matrix M to the type of matrix P and put the result in A:A=typeofMatrix(P)(M)"
},

{
    "location": "MainModule/#Threads-1",
    "page": "MainModule (PosDefManifold.jl)",
    "title": "Threads",
    "category": "section",
    "text": "Some functions in PosDefManifold explicitly call BLAS routines for optimal performnce. This is reported in the help section of the concerned functions. Most functions calls BLAS routine implicitly via Julia. You can set the number of threads the BLAS library should use by:using LinearAlgebra\nBLAS.set_num_threads(n)where n is the number of threads. By default, PosDefManifold reserves to BLAS all CPU threads available on your computer (given by the output of Sys.CPU_THREADS) minus the number threads used by Julia for multi-threaded computations (given by the output of Threads.nthreads()). In Windows this latter number of threads is set to hals the avaiable threads. In Linux and OSX defaults to one and is controlled by an environment variable, i.e.,export JULIA_NUM_THREADS=4See this post, this post and julia doc on threads.Notice that PosDefManifold features many multi-threaded functions and these may allow a gain in computation time only if Julia is instructed to use at least two threads."
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
    "location": "riemannianGeometry/#PosDefManifold.geodesic",
    "page": "riemannianGeometry.jl",
    "title": "PosDefManifold.geodesic",
    "category": "function",
    "text": "(1) geodesic(metric::Metric, P::ℍ{T}, Q::ℍ{T}, a::Real) where T<:RealOrComplex\n(2) geodesic(metric::Metric, D::𝔻{S}, E::𝔻{S}, a::Real) where S<:Real\n\n(1) Move along the geodesic from point P to point Q  (two positive definite matrices) with arclegth 0=a=1,  using the specified metric, of type Metric::Enumerated type.\n\nFor all metrics,\n\nwith a=0 we stay at P,\nwith a=1 we move up to Q,\nwith a=12 we move to the mid-point of P and Q (mean).\n\nUsing the Fisher metric, argument a can be any real number, for instance:\n\nwith 0a1 we move toward Q (attraction),\nwith a1 we move over and beyond Q (extrapolation),\nwith a0 we move back away from Q (repulsion).\n\nP and Q must be flagged by julia as Hermitian.  See typecasting matrices.\n\nNote that if Q=I, the Fisher geodesic move is simply P^a  (no need to call this funtion then).\n\nnote: Nota Bene\nFor the logdet zero and Jeffrey metric no closed form expression for the geodesic is available to the best of authors\' knowledge, so in this case the geodesic is found as the weighted mean using mean. For the Von Neumann not even an expression for the mean is available, so in this case the geodesic is not provided and a warning is printed.\n\n(2) Like in (1), but for two real positive definite diagonal matrices  D and E.\n\nMaths\n\nFor points P, Q and arclength a, letting b=1-a,  the geodesic equations for the supported metrics are:\n\nMetric geodesic equation\nEuclidean bP + aQ\ninvEuclidean big(bP^-1 + aQ^-1big)^-1\nChoEuclidean TT^*, where T=bL_P + aL_Q\nlogEuclidean textexpbig(bhspace2pttextlog(P) + ahspace2pttextlog(Q)big)\nlogCholesky TT^*, where T=S_P+a(S_Q-S_P)+D_Phspace2pttextexpbig(a(textlogD_Q-textlogD_P)big)\nFisher P^12 big(P^-12 Q P^-12big)^a P^12\nlogdet0 uses weighted mean algorithm logdet0Mean\nJeffrey uses weighted mean mean\nVonNeumann N.A.\nWasserstein b^2P+a^2Q +abbig(PQ)^12 +(QP)^12big\n\nlegend: L_X, S_X and D_X    are the Cholesky lower triangle of X, its strictly lower triangular part    and diagonal part, respectively (hence, S_X+D_X=L_X,  L_XL_X^*=X).\n\nSee also: mean.\n\nExamples\n\nusing PosDefManifold\nP=randP(10)\nQ=randP(10)\n# Wasserstein mean\nM=geodesic(Wasserstein, P, Q, 0.5)\n# extrapolate suing the Fisher metric\nE=geodesic(Fisher, P, Q, 2)\n\n\n\n\n\n"
},

{
    "location": "riemannianGeometry/#Geodesic-equations-1",
    "page": "riemannianGeometry.jl",
    "title": "Geodesic equations",
    "category": "section",
    "text": "Function Description\ngeodesic Geodesic equations (weighted mean of two positive definite matrices) for any metric⋅geodesic"
},

{
    "location": "riemannianGeometry/#PosDefManifold.distanceSqr",
    "page": "riemannianGeometry.jl",
    "title": "PosDefManifold.distanceSqr",
    "category": "function",
    "text": "(1) distanceSqr(metric::Metric, P::ℍ{T}) where T<:RealOrComplex\n(2) distanceSqr(metric::Metric, P::ℍ{T}, Q::ℍ{T}) where T<:RealOrComplex\n(3) distanceSqr(metric::Metric, D::𝔻{S}) where S<:Real\n(4) distanceSqr(metric::Metric, D::𝔻{S}, E::𝔻{S}) where S<:Real\n\nalias: distance²\n\n(1) Return δ^2(P I), the square of the distance (or divergence) of positive definite  matrix P from the the identity matrix. See distance from the origin.\n\n(2) Return δ^2(P Q), the square of the distance (or divergence) between two  positive definite matrices P and Q. See distance.\n\nIn both cases the distance function δ is induced by the argument metric of type  Metric::Enumerated type.\n\nP in (1) and P, Q in (2) must be flagged by julia as Hermitian.  See typecasting matrices.\n\n(3) and (4) are specialized methods of (1) and (2), respectively,  for real positive definite Diagonal matrices.  See ℍVector type and 𝔻Vector type.\n\nMaths\n\nFor point P the squared distances from the identity  for the supported metrics are:\n\nMetric Squared Distance from the identity\nEuclidean P-I^2\ninvEuclidean P^-1-I^2\nChoEuclidean L_P-I^2\nlogEuclidean textrmlogP^2\nlogCholesky S_P^2+textrmlogD_P^2\nFisher textrmlogP^2\nlogdet0 textrmlogdetfrac12(P+I) - frac12textrmlogdet(P)\nJeffrey frac12textrmtr(P+P^-1)-n\nVonNeumann frac12textrmtr(PtextrmlogP-textrmlogP)\nWasserstein textrmtr(P+I) -2textrmtr(P^12)\n\nFor points P and Q their squared distances for the supported metrics are:\n\nMetric Squared Distance\nEuclidean P-Q^2\ninvEuclidean P^-1-Q^-1^2\nChoEuclidean  L_P - L_Q ^2\nlogEuclidean textrmlogP-textrmlogQ^2\nlogCholesky S_P-S_Q^2+textrmlogD_P-textrmlogD_Q^2\nFisher textrmlog(P^-12QP^-12)^2\nlogdet0 textrmlogdetfrac12(P+Q) - frac12textrmlogdet(PQ)\nJeffrey frac12textrmtr(Q^-1P+P^-1Q)-n\nVonNeumann frac12textrmtr(PtextrmlogP-PtextrmlogQ+QtextrmlogQ-QtextrmlogP)\nWasserstein textrmtr(P+Q) -2textrmtr(P^12QP^12)^12\n\nlegend: L_X, S_X and D_X   are the Cholesky lower triangle of X, its strictly lower triangular part   and diagonal part, respectively (hence, S_X+D_X=L_X,  L_XL_X^*=X).\n\nSee also: distanceSqrMat.\n\nExamples (1)\n\nusing PosDefManifold\nP=randP(10)\nd=distanceSqr(Wasserstein, P)\ne=distanceSqr(Fisher, P)\nmetric=Metric(Int(logdet0)) # or metric=logdet0\ns=string(metric) # check what is the current metric\nf=distance²(metric, P) #using the alias distance²\n\nExamples (2)\n\nusing PosDefManifold\nP=randP(10)\nQ=randP(10)\nd=distanceSqr(logEuclidean, P, Q)\ne=distance²(Jeffrey, P, Q)\n\n\n\n\n\n"
},

{
    "location": "riemannianGeometry/#PosDefManifold.distance",
    "page": "riemannianGeometry.jl",
    "title": "PosDefManifold.distance",
    "category": "function",
    "text": "(1) distance(metric::Metric, P::ℍ{T}) where T<:RealOrComplex\n(2) distance(metric::Metric, P::ℍ{T}, Q::ℍ{T}) where T<:RealOrComplex\n(3) distance(metric::Metric, D::𝔻{S}) where S<:Real\n(4) distance(metric::Metric, D::𝔻{S}, E::𝔻{S}) where S<:Real\n\n(1) Return δ(P I), the distance between positive definite matrix P and  the identity matrix.\n\n(2) Return δ(P Q), the distance between positive definite  matrices P and Q.\n\n(3) and (4) are specialized methods of (1) and (2), respectively,  for real positive definite Diagonal matrices.\n\nThis is the square root of distanceSqr  and is invoked with the same syntax therein.\n\nSee also: distanceMat.\n\n\n\n\n\n"
},

{
    "location": "riemannianGeometry/#Distances-1",
    "page": "riemannianGeometry.jl",
    "title": "Distances",
    "category": "section",
    "text": "Function Description\ndistanceSqr, distance² Squared distance between positive definite matrices\ndistance Distance between positive definite matrices⋅distanceSqr\ndistance"
},

{
    "location": "riemannianGeometry/#PosDefManifold.distanceSqrMat",
    "page": "riemannianGeometry.jl",
    "title": "PosDefManifold.distanceSqrMat",
    "category": "function",
    "text": "(1) distanceSqrMat(metric::Metric, 𝐏::ℍVector;\n                                    <⏩=false>)\n(2) distanceSqrMat(type::Type{T}, metric::Metric, 𝐏::ℍVector;\n                                    <⏩=false>) where T<:AbstractFloat\n\nalias: distance²Mat\n\nGiven a 1d array 𝐏 of k positive definite matrices  P_1P_k of ℍVector type, create the kk real  LowerTriangular matrix comprising elements δ^2(P_i P_j)textrm for all i=j.\n\nThis is the lower triangular matrix holding all squared inter-distances  (zero on diagonal), using the  specified metric, of type Metric::Enumerated type,  giving rise to distance function δ. See distanceSqr.\n\nOnly the lower triangular part is computed in order to optimize memory use.\n\nBy default, the result matrix is of type Float32. The type can be changed  to another real type using method (2).\n\n<optional keyword arguments>:\n\nif ⏩=true the computation of inter-distances is multi-threaded.\n\nwarning: Multi-Threading\nMulti-threading is still experimental in julia. You should check the result on each computer. Multi-threading is automatically disabled if the number of threads Julia is instructed to use is 2 or 2k. See Threads.\n\nSee: distance.\n\nSee also: laplacian, laplacianEigenMaps, spectralEmbedding.\n\nExamples\n\nusing PosDefManifold\n# Generate a set of 8 random 10x10 SPD matrices\nPset=randP(10, 8) # or, using unicode: 𝐏=randP(10, 8)\n# Compute the squared inter-distance matrix according to the log Euclidean metric.\n# This is much faster as compared to the Fisher metric and in general\n# it is a good approximation.\nDsqr=distanceSqrMat(logEuclidean, Pset)\n# or, using unicode: Δ²=distanceSqrMat(logEuclidean, 𝐏)\n\n# return a matrix of type Float64\nDsqr64=distanceSqrMat(Float64, logEuclidean, Pset)\n\n# Multi-threaded\nDsqr=distanceSqrMat(Fisher, Pset; ⏩=true)\n\n\n\n\n\n"
},

{
    "location": "riemannianGeometry/#PosDefManifold.distanceMat",
    "page": "riemannianGeometry.jl",
    "title": "PosDefManifold.distanceMat",
    "category": "function",
    "text": "(1) distanceMat(metric::Metric, 𝐏::ℍVector;\n                                <⏩=true>)\n(2) distanceMat(type::Type{T}, metric::Metric, 𝐏::ℍVector;\n                                <⏩=true>) where T<:AbstractFloat\n\nGiven a 1d array 𝐏 of k positive definite matrices  P_1P_k of ℍVector type, create the kk real  LowerTriangular matrix comprising elements  δ(P_i P_j)textrm for all i=j.\n\nThis is the lower triangular matrix holding all inter-distances  (zero on diagonal), using the  specified metric, of type Metric::Enumerated type,  giving rise to distance δ. See distance.\n\nOnly the lower triangular part is computed in order to optimize memory use.\n\nBy default, the result matrix is of type Float32. The type can be changed  to another real type using method (2).\n\nThe elements of this matrix are the square root of  distanceSqrMat.\n\n<optional keyword arguments>:\n\nif ⏩=true the computation of inter-distances is multi-threaded.\n\nwarning: Multi-Threading\nMulti-threading is still experimental in julia. You should check the result on each computer. Multi-threading is automatically disabled if the number of threads Julia is instructed to use is 2 or 4k. See Threads.\n\nSee: distance.\n\nExamples\n\nusing PosDefManifold\n# Generate a set of 4 random 10x10 SPD matrices\nPset=randP(10, 4) # or, using unicode: 𝐏=randP(10, 4)\nD=distanceMat(Fisher, Pset)\n# or, using unicode: Δ=distanceMat(Fisher, 𝐏)\n\n# return a matrix of type Float64\nD64=distanceMat(Float64, Fisher, Pset)\n\n# Multi-threaded\nD64=distanceMat(Fisher, Pset; ⏩=true)\n\n\n\n\n\n"
},

{
    "location": "riemannianGeometry/#PosDefManifold.laplacian",
    "page": "riemannianGeometry.jl",
    "title": "PosDefManifold.laplacian",
    "category": "function",
    "text": "laplacian(Δ²:𝕃{S}) where S<:Real\n\nGiven a LowerTriangular matrix of squared inter-distances Δ^2,  return the lower triangular part of the normalized Laplacian.  The elements of the Laplacian are of the same type as the elements of Δ^2.  The result is a LowerTriangular matrix.\n\nFirst, a Gaussian radial basis functions  is applied to all elements of Δ^2, such as\n\nW_ij = exp(fracdisplaystyle-Δ^2_ijdisplaystyleε),\n\nwhere ε is the Gaussian scale parameter chosen automatically   as the median of the elements Δ^2_ij.\n\nFinally, the normalized Laplacian is defined as\n\nΩ = D^-12WD^-12,\n\nwhere D is the diagonal matrix holding on the main diagonal   the sum of the rows (or columns) of W.\n\nnote: Nota Bene\nThe normalized Laplacian as here defined can be requested for any input matrix of squared inter-distances, for example, those obtained on scalars or on vectors using appropriate metrics. In any case, only the lower triangular part of the Laplacian is taken as input. See typecasting matrices.\n\nSee also: distanceSqrMat, laplacianEigenMaps, spectralEmbedding.\n\nExamples\n\nusing PosDefManifold\n# Generate a set of 4 random 10x10 SPD matrices\nPset=randP(10, 4) # or, using unicode: 𝐏=randP(10, 4)\nDsqr=distanceSqrMat(Fisher, Pset) # or: Δ²=distanceSqrMat(Fisher, 𝐏)\nlap=laplacian(Dsqr) # or: Ω=laplacian(Δ²)\n\n\n\n\n\n"
},

{
    "location": "riemannianGeometry/#PosDefManifold.laplacianEigenMaps",
    "page": "riemannianGeometry.jl",
    "title": "PosDefManifold.laplacianEigenMaps",
    "category": "function",
    "text": "laplacianEigenMaps(Ω::𝕃{S}, q::Int;\n                  <tol::Real=0, maxiter=300, ⍰=false>) where S<:Real\n\nalias: laplacianEM\n\nGiven the lower triangular part of a normalized Laplacian Ω  (see laplacian ) return the eigen maps in q dimensions,  i.e., the q eigenvectors of  the normalized Laplacian associated with the largest q  eigenvalues, excluding the first (which is always equal to 1.0).  The eigenvectors are of the same type as Ω.\n\nThe eigenvectors of the normalized Laplacian are computed by the  power iterations+modified Gram-Schmidt method,  allowing the execution of this function for big Laplacian matrices.\n\nReturn the 4-tuple (Λ U iterations convergence), where:\n\nΛ is a qq diagonal matrix holding on diagonal the eigenvalues corresponding to the q dimensions of the Laplacian eigen maps,\nU holds in columns the eigen maps, that is, the q eigenvectors,\niterations is the number of iterations executed by the power method,\nconvergence is the convergence attained by the power method.\n\nThe eigenvectors of U holds the coordinates of the points in a  low-dimension Euclidean space (typically two or three).  This is done for, among other purposes, classifying them and  following their trajectories over time or other dimensions.  For examples of applications see Ridrigues et al. (2018) 🎓  and references therein.\n\nArguments: (Ω::𝕃, q; <tol::Real=0, maxiter=300, ⍰=false>):\n\nΩ is a real LowerTriangular normalized Laplacian obtained by the laplacian function,\nq is the dimension of the Laplacian eigen maps;\nThe following are <optional keyword arguments> for the power iterations:\ntol is the tolerance for convergence (see below),\nmaxiter is the maximum number of iterations allowed,\nif ⍰ is true, the convergence at all iterations will be printed.\n\nnote: Nota Bene\nThe maximum value of q that can be requested is n-1, where n is the size of the Laplacian. In general, q=2 or q=3 is requested.tol defaults to the square root of Base.eps of the (real) type of Ω. This corresponds to requiring equality for the convergence criterion over two successive power iterations of about half of the significant digits.\n\nSee also: distanceSqrMat, laplacian, spectralEmbedding.\n\nExamples\n\nusing PosDefManifold\n# Generate a set of 4 random 10x10 SPD matrices\nPset=randP(10, 4) # or, using unicode: 𝐏=randP(10, 4)\nDsqr=distanceSqrMat(Fisher, Pset) #or: Δ²=distanceSqrMat(Fisher, 𝐏)\nlap= laplacian(Dsqr) # or: Ω=laplacian(Δ²)\nevalues, maps, iterations, convergence=laplacianEM(lap, 2)\nevalues, maps, iterations, convergence=laplacianEM(lap, 2; maxiter=100)\nevalues, maps, iterations, convergence=laplacianEM(lap, 2; ⍰=true)\n\n\n\n\n\n"
},

{
    "location": "riemannianGeometry/#PosDefManifold.spectralEmbedding",
    "page": "riemannianGeometry.jl",
    "title": "PosDefManifold.spectralEmbedding",
    "category": "function",
    "text": "(1) spectralEmbedding(metric::Metric, 𝐏::ℍVector, q::Int;\n                     <tol::Real=0, maxiter=300, ⍰=false>)\n\n(2) spectralEmbedding(metric::Metric, 𝐏::ℍVector, q::Int, type::Type{T};\n                     <tol::Real=0, maxiter=300, ⍰=false>) where T<:Real\n\nGiven a 1d array 𝐏 of k positive definite matrices P_1P_k  (real or complex), compute its eigen maps in q dimensions.\n\nThis function runs one after the other the functions:\n\ndistanceSqrMat (compute the squared inter-distance matrix),\nlaplacian (compute the normalized Laplacian),\nlaplacianEigenMaps (get the eigen maps).\n\nBy default all computations above are done with Float32 precision.  Another real type can be requested using method (2), where the type argument  is defined.\n\nReturn the 4-tuple (Λ, U, iterations, convergence), where:\n\nΛ is a qq diagonal matrix holding on diagonal the eigenvalues corresponding to the q dimensions of the Laplacian eigen maps,\nU holds in columns the q eigenvectors, i.e., the q coordinates of the points in the embedded space,\niterations is the number of iterations executed by the power method,\nconvergence is the convergence attained by the power method.\n\nArguments (metric, 𝐏, q, <tol::Real=0, maxiter=300, ⍰=false>):\n\nmetric is the metric of type Metric::Enumerated type used for computing the inter-distances,\n𝐏 is a 1d array of k positive matrices of ℍVector type,\nq is the dimension of the Laplacian eigen maps;\nThe following are <optional keyword arguments> for the power method iterative algorithm:\ntol is the tolerance for convergence of the power method (see below),\nmaxiter is the maximum number of iterations allowed for the power method,\nif ⍰ is true the convergence at all iterations will be printed.\n\nnote: Nota Bene\ntol defaults to the square root of Base.eps of the Float32 type (1) or of the type passed as argumant (2). This corresponds to requiring equality for the convergence criterion over two successive power iterations of about half of the significant digits.\n\nSee also: distanceSqrMat, laplacian, laplacianEigenMaps.\n\nExamples\n\nusing PosDefManifold\n# Generate a set of 4 random 10x10 SPD matrices\nPset=randP(10, 4) # or, using unicode: 𝐏=randP(10, 4)\nevalues, maps, iterations, convergence=spectralEmbedding(logEuclidean, Pset, 2)\n# show convergence information\nevalues, maps, iterations, convergence=spectralEmbedding(logEuclidean, Pset, 2; ⍰=true)\n# use Float64 precision.\nevalues, maps, iterations, convergence=spectralEmbedding(logEuclidean, Pset, 2, Float64)\n\n\n\n\n\n"
},

{
    "location": "riemannianGeometry/#Graphs-and-Laplacians-1",
    "page": "riemannianGeometry.jl",
    "title": "Graphs and Laplacians",
    "category": "section",
    "text": "Function Description\ndistanceSqrMat, distance²Mat Lower triangular matrix of all squared inter-distances\ndistanceMat Lower triangular matrix of all inter-distances\nlaplacian Laplacian of a squared inter-distances matrix\nlaplacianEigenMaps, laplacianEM Eigen maps (eigenvectors) of a Laplacian\nspectralEmbedding Spectral Embedding (all the above function run in series)⋅distanceSqrMat\ndistanceMat\nlaplacian\nlaplacianEigenMaps\nspectralEmbedding"
},

{
    "location": "riemannianGeometry/#Statistics.mean",
    "page": "riemannianGeometry.jl",
    "title": "Statistics.mean",
    "category": "function",
    "text": "(1) mean(metric::Metric, P::ℍ{T}, Q::ℍ{T}) where T<:RealOrComplex\n(2) mean(metric::Metric, D::𝔻{T}, E::𝔻{T}) where T<:Real\n\n(3) mean(metric::Metric, 𝐏::ℍVector;\n        <w::Vector=[], ✓w=true, ⏩=false>)\n(4) mean(metric::Metric, 𝐃::𝔻Vector;\n        <w::Vector=[], ✓w=true, ⏩=false>)\n\n(1) Mean of two positive definite matrices, passed in arbitrary order as  arguments P and Q, using the specified metric of type  Metric::Enumerated type.  The order is arbitrary as all metrics implemented in PosDefManifold are symmetric.  This is the midpoint of the geodesic.  For the weighted mean of two positive definite matrices use instead  the geodesic function.  P and Q must be flagged as Hermitian. See typecasting matrices.\n\n(2) Like in (1), but for two real diagonal positive definite matrices  D and E.\n\n(3) Fréchet mean of an 1d array 𝐏 of k positive definite  matrices 𝐏=P_1P_k of ℍVector type,  with optional non-negative real weights w=w_1w_k and using the  specified metricas in (1).\n\n(4) Fréchet mean of an 1d array 𝐃 of k positive definite  matrices 𝐃=D_1D_k of 𝔻Vector type,  with optional non-negative real weights w=w_1w_k and using the  specified metricas in (1).\n\nIf you don\'t pass a weight vector with <optional keyword argument> w,  return the unweighted mean.\n\nIf <optional keyword argument> ✓w=true (default), the weights are  normalized so as to sum up to 1, otherwise they are used as they are passed  and should be already normalized.  This option is provided to allow  calling this function repeatedly without normalizing the same weights  vector each time.\n\nAdopting the Fisher, logdet0 and Wassersteinmetric in (3) and thelogdet0metric in (4), the mean is computed by means of an iterative  algorithm and information on its convergence is displayed in the REPL.  For suppressing this information and for more options for computing these means  call directly functions [geometricMean](@ref), [logdet0Mean](@ref)  and [wasMean`](@ref).\n\nfor (3) and (4), if ⏩=true the computation of the mean is multi-threaded.\n\nwarning: Multi-Threading\nMulti-threading is still experimental in julia. You should check the result on each computer. Multi-threading is automatically disabled if the number of threads Julia is instructed to use is 2 or 4k. See Threads.\n\nMath\n\nThe Fréchet mean of a set of k matrices P_1 P_2 P_k weighted by  w_1 w_2 w_ksum_i=1^kw_i=1 for the supported metrics are,  for those with closed form expression:\n\nMetric weighted Fréchet mean\nEuclidean sum_i=1^kw_i P_i\ninvEuclidean big(sum_i=1^kw_i P_i^-1big)^-1\nChoEuclidean TT^*, where T=bL_P + aL_Q\nlogEuclidean textrmexpbig(sum_i=1^kw_ihspace1pt textrmlogP_i big)\nlogCholesky TT^*, where T=sum_i=1^k(w_kS_k)+sum_i=1^k(w_ktextrmlogD_k)\nJeffrey A^12big(A^-12HA^-12big)^12A^12\n\nand for those that verify an equation:\n\nMetric equation verified by the weighted Fréchet mean\nFisher sum_i=1^kw_itextrmlogbig(G^-12 P_k G^-12big)=0\nlogdet0 sum_i=1^kw_ibig(frac12P_i+frac12Gbig)^-1=G^-1\nVonNeumann N.A.\nWasserstein G=sum_i=1^kw_ibig( G^12  P_i G^12big)^12\n\nlegend: L_X, S_X and D_X   are the Cholesky lower triangle of X, its strictly lower triangular part   and diagonal part, respectively (hence, S_X+D_X=L_X,  L_XL_X^*=X).   A and H are the weighted arithmetic and weighted harmonic mean, respectively.\n\nSee: geodesic, mean, Fréchet mean.\n\nExamples\n\nusing LinearAlgebra, Statistics, PosDefManifold\n# Generate 2 random 3x3 SPD matrices\nP=randP(3)\nQ=randP(3)\nM=mean(logdet0, P, Q) # (1)\nM=mean(logdet0, P, Q) # (1)\n\nR=randP(3)\n# passing several matrices and associated weights listing them\n# weights vector, does not need to be normalized\nmean(Fisher, ℍVector([P, Q, R]); w=[1, 2, 3])\n\n# Generate a set of 4 random 3x3 SPD matrices\nPset=randP(3, 4) # or, using unicode: 𝐏=randP(3, 4)\nweights=[1, 2, 3, 1]\n# passing a vector of Hermitian matrices (ℍVector type)\nM=mean(Euclidean, Pset; w=weights) # (2) weighted Euclidean mean\nM=mean(Wasserstein, Pset)  # (2) unweighted Wassertein mean\n# using unicode: M=mean(Wasserstein, 𝐏)\n\n# run multi-threaded when the number of matrices is high\nusing BenchmarkTools\nPset=randP(20, 160)\n@benchmark(mean(logEuclidean, Pset)) # single-threaded\n@benchmark(mean(logEuclidean, Pset; ⏩=true)) # multi-threaded\n\n\n\n\n\n"
},

{
    "location": "riemannianGeometry/#PosDefManifold.means",
    "page": "riemannianGeometry.jl",
    "title": "PosDefManifold.means",
    "category": "function",
    "text": "(1) means(metric::Metric, 𝒫::ℍVector₂; <⏩=false>)\n\n(2) means(metric::Metric, 𝒟::𝔻Vector₂; <⏩=false>)\n\n(1) Given a 2d array 𝒫 of positive definite matrices as an ℍVector₂ type  compute the Fréchet mean for as many ℍVector type objects  as hold in 𝒫, using the specified metric of type  Metric::Enumerated type.  Return the means in a vector of Hermitian matrices, that is, as an ℍVector type.\n\n(2) Given a 2d array 𝒟 of real positive definite matrices as an 𝔻Vector₂ type  compute the Fréchet mean for as many 𝔻Vector type objects  as hold in 𝒟, using the specified metric of type  Metric::Enumerated type.  Return the means in a vector of Diagonal matrices, that is, as a 𝔻Vector type.\n\nThe weigted Fréchet mean is not supported in this function.\n\nIf <optional key argmuent> ⏩=true the computation of the means  is multi-threaded.\n\nwarning: Multi-Threading\nMulti-threading is still experimental in julia. You should check the result on each computer. For each mean to be computed, multi-threading is automatically disabled if the number of threads Julia is instructed to use is 2 or 4k, where k is the number of matrices for which the mean is to be computed. See Threads.\n\nSee also: mean.\n\nExamples\n\n using PosDefManifold\n # Generate a set of 4 random 3x3 SPD matrices\n Pset=randP(3, 4) # or, using unicode: 𝐏=randP(3, 4)\n # Generate a set of 40 random 4x4 SPD matrices\n Qset=randP(3, 40) # or, using unicode: 𝐐=randP(3, 40)\n # listing directly ℍVector objects\n means(logEuclidean, ℍVector₂([Pset, Qset])) # or: means(logEuclidean, ℍVector₂([𝐏, 𝐐]))\n # note that [𝐏, 𝐐] is actually a ℍVector₂ type object\n\n # creating and passing an object of ℍVector₂ type\n sets=ℍVector₂(undef, 2) # or: 𝒫=ℍVector₂(undef, 2)\n sets[1]=Pset # or: 𝒫[1]=𝐏\n sets[2]=Qset # or: 𝒫[2]=𝐐\n means(logEuclidean, sets) # or: means(logEuclidean, 𝒫)\n\n # going multi-threated\n\n # first, create 20 sets of 200 50x50 SPD matrices\n sets=ℍVector₂([randP(50, 200) for i=1:20])\n\n # How much computing time we save ?\n # (example min time obtained with 4 threads & 4 BLAS threads)\n using BenchmarkTools\n\n # non multi-threaded, mean with closed-form solution\n @benchmark(means(logEuclidean, sets))  		 # (6.196 s)\n\n # multi-threaded, mean with closed-form solution\n @benchmark(means(logEuclidean, sets; ⏩=true)) # (1.897 s)\n\n sets=ℍVector₂([randP(10, 200) for i=1:10])\n\n # non multi-threaded, mean with iterative solution\n # wait a bit\n @benchmark(means(Fisher, sets))  		         # (4.672 s )\n\n # multi-threaded, mean with iterative solution\n @benchmark(means(Fisher, sets; ⏩=true))        # (1.510 s)\n\n\n\n\n\n"
},

{
    "location": "riemannianGeometry/#PosDefManifold.generalizedMean",
    "page": "riemannianGeometry.jl",
    "title": "PosDefManifold.generalizedMean",
    "category": "function",
    "text": "generalizedMean(𝐏::Union{ℍVector, 𝔻Vector}, p::Real;\n               <w::Vector=[], ✓w=true, ⏩=false>)\n\nGiven a 1d array 𝐏=P_1P_k of k positive definite matrices of  ℍVector type or real positive definite diagonal matrices of  𝔻Vector type and optional non-negative real weights vector  w=w_1w_k, return the weighted generalized means G  with real parameter p, that is,\n\nG=big(sum_i=1^kw_iP_i^pbig)^1p.\n\nIf you don\'t pass a weight vector with <optional keyword argument> w,  return the unweighted generalized mean\n\nG=big(sum_i=1^kP_i^pbig)^1p.\n\nIf <optional keyword argument> ✓w=true (default), the weights are  normalized so as to sum up to 1, otherwise they are used as they are passed  and should be already normalized.  This option is provided to allow  calling this function repeatedly without normalizing the weights each time.\n\nIf <optional key argmuent> ⏩=true the computation of the generalized mean  is multi-threaded.\n\nwarning: Multi-Threading\nMulti-threading is still experimental in julia. You should check the result on each computer. Multi-threading is automatically disabled if the number of threads Julia is instructed to use is 2 or 4k. See Threads.\n\nThe following special cases for parameter p are noteworthy:\n\nFor p=frac12 the generalized mean is the modified Bhattacharyya mean.\nFor p=1 the generalized mean is the Euclidean mean.\nFor p=-1 the generalized mean is the inverse Euclidean mean.\nFor p=0 the generalized mean is the log Euclidean mean, which is the Fisher mean when matrices in 𝐏 all pair-wise commute.\n\nNotice that when matrices in 𝐏 all pair-wise commute, for instance if the  matrices are diagonal,  the generalized means coincide with the power means  for any p-1 1 and for p=05 it coincides also with the  Wasserstein mean. For this reason the generalized means are used  as default initialization of both the powerMean and wasMean  algorithm.\n\nSee: generalized means.\n\nSee also: powerMean, wasMean, mean.\n\nExamples\n\nusing LinearAlgebra, Statistics, PosDefManifold\n# Generate a set of 4 random 3x3 SPD matrices\nPset=randP(3, 4) # or, using unicode: 𝐏=randP(3, 4)\n\n# weights vector, does not need to be normalized\nweights=[1, 2, 3, 1]\n\n# unweighted mean\nG = generalizedMean(Pset, 0.25) # or: G = generalizedMean(𝐏, 0.25)\n\n# weighted mean\nG = generalizedMean(Pset, 0.5; w=weights)\n\n# with weights previously normalized we can set ✓w=false\nweights=weights./sum(weights)\nG = generalizedMean(Pset, 0.5; w=weights, ✓w=false)\n\n# run multi-threaded when the number of matrices is high\nusing BenchmarkTools\nPset=randP(20, 160)\n@benchmark(generalizedMean(Pset)) # single-threaded\n@benchmark(generalizedMean(Pset; ⏩=true)) # multi-threaded\n\n\n\n\n\n"
},

{
    "location": "riemannianGeometry/#PosDefManifold.geometricMean",
    "page": "riemannianGeometry.jl",
    "title": "PosDefManifold.geometricMean",
    "category": "function",
    "text": "geometricMean(𝐏::Union{ℍVector, 𝔻Vector};\n        <w::Vector=[], ✓w=true, init=nothing, tol::Real=0, ⍰=false, ⏩=false>)\n\nalias: gmean\n\nGiven a 1d array 𝐏=P_1P_k of k positive definite matrices of  ℍVector type or diagonal matrices of 𝔻Vector type  and optional non-negative real weights vector w=w_1w_k,  return the 3-tuple (G iter conv), where G is the mean according  to the Fisher metric and iter, conv are the number of iterations  and convergence attained by the algorithm.  Mean G is the unique positive definite matrix satisfying\n\nsum_i=1^kw_itextrmlogbig(G^-12 P_i G^-12big)=0\n\nFor estimating it, this function implements the well-known gradient descent  algorithm, yielding iterations\n\nG G^12textrmexpbig(sum_i=1^kw_itextrmlog(G^-12 P_i G^-12)big)G^12\n\nIf you don\'t pass a weight vector with <optional keyword argument> w,  return the unweighted geometric mean.\n\nIf <optional keyword argument> ✓w=true (default), the weights are  normalized so as to sum up to 1, otherwise they are used as they are passed  and should be already normalized.  This option is provided to allow  calling this function repeatedly without normalizing the same weights  vector each time.\n\nThe following are more <optional keyword arguments>:\n\ninit is a matrix to be used as initialization for the mean. If no matrix is provided, the log Euclidean mean will be used,\ntol is the tolerance for the convergence (see below).\nif ⍰ is true, the convergence attained at each iteration is printed.\nif ⏩=true the iterations are multi-threaded.\n\nIf the input is a 1d array of k real positive definite diagonal matrices  the solution is available in closed-form as the log Euclidean  mean, hence the <optional keyword arguments> init, tol and ⍰  have no effect and return the 3-tuple (G 1 0).  See the log Euclidean metric.\n\nwarning: Multi-Threading\nMulti-threading is still experimental in julia. You should check the result on each computer. Multi-threading is automatically disabled if the number of threads Julia is instructed to use is 2 or 4k. See Threads.\n\nnote: Nota Bene\nIn normal circumstances this algorithm converges monothonically. If the algorithm diverges a warning is printed indicating the iteration when this happened.tol defaults to 100 times the square root of Base.eps of the nearest real type of data input 𝐏. This corresponds to requiring the relative convergence criterion over two successive iterations to vanish for about half the significant digits minus 2.\n\nSee: Fisher metric.\n\nSee also: powerMean, wasMean, logdet0Mean,  mean.\n\nExamples\n\nusing LinearAlgebra, PosDefManifold\n# Generate a set of 4 random 3x3 SPD matrices\nPset=randP(3, 4) # or, using unicode: 𝐏=randP(3, 4)\n\n# unweighted mean\nG, iter, conv = geometricMean(Pset) # or G, iter, conv = geometricMean(𝐏)\n\n# weights vector, does not need to be normalized\nweights=[1, 2, 3, 1]\n\n# weighted mean\nG, iter, conv = geometricMean(Pset, w=weights)\n\n# print the convergence at all iterations\nG, iter, conv = geometricMean(Pset; w=weights, ⍰=true)\n\n# now suppose Pset has changed a bit, initialize with G to hasten convergence\nPset[1]=ℍ(Pset[1]+(randP(3)/100))\nG, iter, conv = geometricMean(Pset; w=weights, ✓w=false, ⍰=true, init=G)\n\n# run multi-threaded when the number of matrices is high\nusing BenchmarkTools\nPset=randP(20, 160)\n@benchmark(geometricMean(Pset)) # single-threaded\n@benchmark(geometricMean(Pset; ⏩=true)) # multi-threaded\n\n\n\n\n\n"
},

{
    "location": "riemannianGeometry/#PosDefManifold.logdet0Mean",
    "page": "riemannianGeometry.jl",
    "title": "PosDefManifold.logdet0Mean",
    "category": "function",
    "text": "logdet0Mean(𝐏::Union{ℍVector, 𝔻Vector};\n    <w::Vector=[], ✓w=true, init=nothing, tol::Real=0, ⍰=false, ⏩=false>)\n\nalias: ld0Mean\n\nGiven a 1d array 𝐏=P_1P_k of k positive definite matrices of  ℍVector type or real positive definite diagonal matrices of  𝔻Vector type and optional  non-negative real weights vector w=w_1w_k,  return the 3-tuple (G iter conv), where G is the mean according  to the logdet zero metric and iter, conv are the number of iterations  and convergence attained by the algorithm.  Mean G is the unique positive definite matrix satisfying\n\nsum_i=1^kw_ibig(frac12P_i+frac12Gbig)^-1=G^-1.\n\nFor estimating it, this function implements the fixed-point iteration algorithm  suggested by (Moakher, 2012, p315)🎓, yielding iterations\n\nG  frac12big(sum_i=1^kw_i(P_i+G)^-1big)^-1.\n\nIf you don\'t pass a weight vector with <optional keyword argument> w,  return the unweighted logdet zero mean.\n\nIf <optional keyword argument> ✓w=true (default), the weights are  normalized so as to sum up to 1, otherwise they are used as they are passed  and should be already normalized.  This option is provided to allow  calling this function repeatedly without normalizing the same weights  vector each time.\n\nThe following are more <optional keyword arguments>:\n\ninit is a matrix to be used as initialization for the mean. If no matrix is provided, the log Euclidean mean will be used,\ntol is the tolerance for the convergence (see below).\nif ⍰ is true, the convergence attained at each iteration is printed.\nif ⏩=true the iterations are multi-threaded.\n\nwarning: Multi-Threading\nMulti-threading is still experimental in julia. You should check the result on each computer. Multi-threading is automatically disabled if the number of threads Julia is instructed to use is 2 or 4k. See Threads.\n\nnote: Nota Bene\nIn normal circumstances this algorithm converges monothonically. If the algorithm diverges a warning is printed indicating the iteration when this happened.tol defaults to 100 times the square root of Base.eps of the nearest real type of data input 𝐏. This corresponds to requiring the relative convergence criterion over two successive iterations to vanish for about half the significant digits minus 2.\n\nSee: logdet zero metric, modified Bhattacharyya mean.\n\nSee also: powerMean, wasMean, logdet0Mean,  mean.\n\nExamples\n\nusing LinearAlgebra, PosDefManifold\n# Generate a set of 4 random 3x3 SPD matrices\nPset=randP(3, 4) # or, using unicode: 𝐏=randP(3, 4)\n\n# unweighted mean\nG, iter, conv = logdet0Mean(Pset) # or G, iter, conv = logdet0Mean(𝐏)\n\n# weights vector, does not need to be normalized\nweights=[1, 2, 3, 1]\n\n# weighted mean\nG, iter, conv = logdet0Mean(Pset, w=weights)\n\n# print the convergence at all iterations\nG, iter, conv = logdet0Mean(Pset; w=weights, ⍰=true)\n\n# now suppose Pset has changed a bit, initialize with G to hasten convergence\nPset[1]=ℍ(Pset[1]+(randP(3)/100))\nG, iter, conv = logdet0Mean(Pset; w=weights, ✓w=false, ⍰=true, init=G)\n\n# run multi-threaded when the number of matrices is high\nusing BenchmarkTools\nPset=randP(20, 160)\n@benchmark(logdet0Mean(Pset)) # single-threaded\n@benchmark(logDet0Mean(Pset; ⏩=true)) # multi-threaded\n\n\n\n\n\n"
},

{
    "location": "riemannianGeometry/#PosDefManifold.wasMean",
    "page": "riemannianGeometry.jl",
    "title": "PosDefManifold.wasMean",
    "category": "function",
    "text": "wasMean(𝐏::Union{ℍVector, 𝔻Vector};\n        <w::Vector=[], ✓w=true, init=nothing, tol::Real=0, ⍰=false, ⏩=false>)\n\nGiven a 1d array 𝐏=P_1P_k of k positive definite matrices  of ℍVector type or real positive definite diagonal matrices of  𝔻Vector type and optional non-negative real weights vector  w=w_1w_k,  return the 3-tuple (G iter conv), where G is the mean according  to the Wasserstein metric and iter, conv are the number of iterations  and convergence attained by the algorithm.  Mean G is the unique positive definite matrix satisfying\n\nG=sum_i=1^kw_ibig( G^12  P_i G^12big)^12.\n\nFor estimating it, this function implements the fixed-point iterative algorithm  proposed by (Álvarez-Esteban et al., 2016)🎓:\n\nG  G^-12big(sum_i=1^k w_i(G^12P_i G^12)^12big)^2 G^-12.\n\nIf you don\'t pass a weight vector with <optional keyword argument> w,  return the unweighted Wassertein mean.\n\nIf <optional keyword argument> ✓w=true (default), the weights are  normalized so as to sum up to 1, otherwise they are used as they are passed  and they should be already normalized.  This option is provided to allow  calling this function repeatedly without normalizing the same weights  vector each time.\n\nThe following are more <optional keyword arguments>:\n\ninit is a matrix to be used as initialization for the mean. If no matrix is provided, the instance of generalized means with p=05 will be used,\ntol is the tolerance for the convergence (see below).\nif ⍰ is true, the convergence attained at each iteration is printed.\nif ⏩=true the iterations are multi-threaded.\n\nIf the input is a 1d array of k real positive definite diagonal matrices  the solution is available in closed-form as the modified Bhattacharyya mean,  hence the <optional keyword arguments> init, tol and ⍰  have no effect and return the 3-tuple (G 1 0).  See modified Bhattacharyya mean.\n\nwarning: Multi-Threading\nMulti-threading is still experimental in julia. You should check the result on each computer. Multi-threading is automatically disabled if the number of threads Julia is instructed to use is 2 or 4k. See Threads.\n\nnote: Nota Bene\nIn normal circumstances this algorithm converges monothonically. If the algorithm diverges a warning is printed indicating the iteration when this happened.tol defaults to 100 times the square root of Base.eps of the nearest real type of data input 𝐏. This corresponds to requiring the relative convergence criterion over two successive iterations to vanish for about half the significant digits minus 2.\n\nSee: Wasserstein metric.\n\nSee also: powerMean, wasMean, logdet0Mean,  mean.\n\nExamples\n\nusing LinearAlgebra, PosDefManifold\n# Generate a set of 4 random 3x3 SPD matrices\nPset=randP(3, 4) # or, using unicode: 𝐏=randP(3, 4)\n\n# unweighted mean\nG, iter, conv = wasMean(Pset) # or: G, iter, conv = wasMean(𝐏)\n\n# weights vector, does not need to be normalized\nweights=[1, 2, 3, 1]\n\n# weighted mean\nG, iter, conv = wasMean(Pset; w=weights)\n\n# print the convergence at all iterations\nG, iter, conv = wasMean(Pset; w=weights, ⍰=true)\n\n# now suppose 𝐏 has changed a bit, initialize with G to hasten convergence\nPset[1]=ℍ(Pset[1]+(randP(3)/100))\nG, iter, conv = wasMean(Pset; w=weights, ⍰=true, init=G)\n\n# run multi-threaded when the number of matrices is high\nusing BenchmarkTools\nPset=randP(20, 160)\n@benchmark(wasMean(Pset)) # single-threaded\n@benchmark(wasMean(Pset; ⏩=true)) # multi-threaded\n\n\n\n\n\n"
},

{
    "location": "riemannianGeometry/#PosDefManifold.powerMean",
    "page": "riemannianGeometry.jl",
    "title": "PosDefManifold.powerMean",
    "category": "function",
    "text": "powerMean(𝐏::Union{ℍVector, 𝔻Vector}, p::Real;\n         <w::Vector=[], ✓w=true, init=nothing, tol::Real=0, ⍰=false, ⏩=false>)\n\nGiven a 1d array 𝐏=P_1P_k of k positive definite matrices  of ℍVector type or real positive definite diagonal matrices of  𝔻Vector type, an optional non-negative real weights vector  w=w_1w_k and a real parameter p in-1 1, return the  3-tuple (G iter conv), where G is  Lim and Palfia (2012)\'s power means  of order p and  iter, conv are the number of iterations  and convergence attained by the algorithm, respectively.  Mean G is the unique positive definite matrix satisfying\n\nG=sum_i=1^k(w_iGtextrm_pP_i),\n\nwhere Gtextrm_pP_i is the Fisher geodesic equation.  In particular:\n\nwith p=-1 this is the harmonic mean (see the inverse Euclidean metric),\nwith p=+1 this is the arithmetic mean (see the Euclidean metric),\nat the limit of p evaluated at zero from both side this is the geometric mean (see Fisher metric).\n\nFor estimating power means for pin(-1 1), this function implements  the  fixed-point iterative algorithm of (Congedo et al., 2017b)🎓.  For p=0 (geometric mean)  this algorithm is run two times with a small positive and negative value  of p and the geometric mean of the two  resulting means is returned, as suggested in (Congedo et al., 2017b)🎓.  This way of estimating the geometric mean of  a set of matrices is faster as compared to the usual gradient descent algorithm.\n\nIf you don\'t pass a weight vector with <optional keyword argument> w,  return the unweighted power mean.\n\nIf <optional keyword argument> ✓w=true (default), the weights are  normalized so as to sum up to 1, otherwise they are used as they are passed  and should be already normalized.  This option is provided to allow  calling this function repeatedly without normalizing the same weights  vector each time.\n\nThe following are more <optional keyword arguments>:\n\ninit is a matrix to be used as initialization for the mean. If no matrix is provided, the instance of generalized means with parameter p will be used.\ntol is the tolerance for the convergence (see below).\nif ⍰ is true, the convergence attained at each iteration is printed.\nif ⏩=true the iterations are multi-threaded.\n\nIf the input is a 1d array of k real positive definite diagonal matrices  the solution is available in closed-form as the generalized  mean of order p, hence the <optional keyword arguments>  init, tol and ⍰  have no effect and return the 3-tuple (G 1 0).  See generalized means.\n\nwarning: Multi-Threading\nMulti-threading is still experimental in julia. You should check the result on each computer. Multi-threading is automatically disabled if the number of threads Julia is instructed to use is 2 or 4k. See Threads.\n\nnote: Nota Bene\nIn normal circumstances this algorithm converges monothonically. If the algorithm diverges a warning is printed indicating the iteration when this happened.tol defaults to 100 times the square root of Base.eps of the nearest real type of data input 𝐏. This corresponds to requiring the relative convergence criterion over two successive iterations to vanish for about half the significant digits minus 2.\n\n(2) Like in (1), but for a 1d array 𝐃=D_1D_k of k  real positive definite diagonal matrices of 𝔻Vector type.  In this case the solution is available in closed-form, hence the  <optional keyword arguments> init, tol and ⍰ have no effect and return  the 3-tuple (G 1 0). See generalized means.\n\nSee: power means, generalized means, modified Bhattacharyya mean.\n\nSee also: generalizedMean, wasMean, logdet0Mean,  mean.\n\nExamples\n\nusing LinearAlgebra, PosDefManifold\n# Generate a set of 4 random 3x3 SPD matrices\nPset=randP(3, 4) # or, using unicode: 𝐏=randP(3, 4)\n\n# unweighted mean\nG, iter, conv = powerMean(Pset, 0.5) # or G, iter, conv = powerMean(𝐏, 0.5)\n\n# weights vector, does not need to be normalized\nweights=[1, 2, 3, 1]\n\n# weighted mean\nG, iter, conv = powerMean(Pset, 0.5; w=weights)\n\n# print the convergence at all iterations\nG, iter, conv = powerMean(Pset, 0.5; w=weights, ⍰=true)\n\n# now suppose 𝐏 has changed a bit, initialize with G to hasten convergence\nPset[1]=ℍ(Pset[1]+(randP(3)/100))\nG, iter, conv = powerMean(Pset, 0.5; w=weights, ⍰=true, init=G)\n\n# run multi-threaded when the number of matrices is high\nusing BenchmarkTools\nPset=randP(20, 160)\n@benchmark(powerMean(Pset, 0.5)) # single-threaded\n@benchmark(powerMean(Pset, 0.5; ⏩=true)) # multi-threaded\n\n\n\n\n\n"
},

{
    "location": "riemannianGeometry/#Means-1",
    "page": "riemannianGeometry.jl",
    "title": "Means",
    "category": "section",
    "text": "Function Description\nmean Weighted Fréchet mean (wFm) of a matrix set using any metric\nmeans As above for several sets at once\ngeneralizedMean Generalized wFm of a matrix set\ngeometricMean, gMean wFm of a matrix set according to the Fisher metric (iterative)\nlogdet0Mean, ld0Mean wFm of a matrix set according to the logdet0 metric (iterative)\nwasMean wFm of a matrix set according to the Wasserstein metric (iterative)\npowerMean Power wFm of a matrix set (iterative)⋅mean\nmeans\ngeneralizedMean\ngeometricMean\nlogdet0Mean\nwasMean\npowerMean"
},

{
    "location": "riemannianGeometry/#PosDefManifold.logMap",
    "page": "riemannianGeometry.jl",
    "title": "PosDefManifold.logMap",
    "category": "function",
    "text": "logMap(metric::Metric, P::ℍ{T}, G::ℍ{T}) where T<:RealOrComplex\n\nLogaritmic Map: map a positive definite matrix P from the SPD or  Hermitian manifold into the tangent space at base-point G using the Fisher metric.\n\nP and G must be flagged as Hermitian. See typecasting matrices.\n\nThe map is defined as\n\nLog_G(P)=S=G^12textrmlogbig(G^-12PG^-12big)G^12.\n\nThe result is an Hermitian matrix.  The inverse operation is expMap.\n\nArguments (metric, P, G):\n\nmetric is a metric of type Metric::Enumerated type.\nP is the positive definite matrix to be projected onto the tangent space,\nG is the tangent space base point,\n\nCurrently only the Fisher metric is supported for tangent space operations.\n\nSee also: vecP.\n\nExamples\n\nusing PosDefManifold\nP=randP(3)\nQ=randP(3)\nmetric=Fisher\nG=mean(metric, P, Q)\n# projecting P at the base point given by the geometric mean of P and Q\nS=logMap(metric, P, G)\n\n\n\n\n\n"
},

{
    "location": "riemannianGeometry/#PosDefManifold.expMap",
    "page": "riemannianGeometry.jl",
    "title": "PosDefManifold.expMap",
    "category": "function",
    "text": "expMap(metric::Metric, S::ℍ{T}, G::ℍ{T}) where T<:RealOrComplex\n\nExponential Map: map an Hermitian matrix S from the tangent space at base  point G into the SPD or Hermitian manifold (using the Fisher metric).\n\nS and G must be flagged as Hermitian. See typecasting matrices.\n\nThe map is defined as\n\nExp_G(S)=P=G^12textrmexpbig(G^-12SG^-12big)G^12.\n\nThe result is a positive definite matrix.  The inverse operation is logMap.\n\nArguments (metric, S, G):\n\nmetric is a metric of type Metric::Enumerated type,\nS is a Hermitian matrix, real or complex, to be projected on the SPD or Hermitian manifold,\nG is the tangent space base point.\n\nCurrently only the Fisher metric is supported for tangent space operations.\n\nExamples\n\nusing PosDefManifold, LinearAlgebra\nP=randP(3)\nQ=randP(3)\nG=mean(Fisher, P, Q)\n# projecting P on the tangent space at the Fisher mean base point G\nS=logMap(Fisher, P, G)\n# adding the identity in the tangent space and reprojecting back onto the manifold\nH=expMap(Fisher, ℍ(S+I), G)\n\n\n\n\n\n"
},

{
    "location": "riemannianGeometry/#PosDefManifold.vecP",
    "page": "riemannianGeometry.jl",
    "title": "PosDefManifold.vecP",
    "category": "function",
    "text": "vecP(S::ℍ{T}) where T<:RealOrComplex\n\nVectorize a tangent vector (which is an Hermitian matrix) S:  mat -> vec.\n\nIt gives weight 1 to diagonal elements and √2 to off-diagonal elements  (Barachant et al., 2012)🎓.\n\nThe result is a vector holding n(n+1)2 elements, where n  is the size of S.\n\nS must be flagged as Hermitian. See typecasting matrices.\n\nThe inverse operation is provided by matP.\n\nExamples\n\nusing PosDefManifold\nP=randP(3)\nQ=randP(3)\nG=mean(Fisher, P, Q)\n# projecting P at the base point given by the geometric mean of P and Q\nS=logMap(Fisher, P, G)\n# vectorize S\nv=vecP(S)\n\n\n\n\n\n"
},

{
    "location": "riemannianGeometry/#PosDefManifold.matP",
    "page": "riemannianGeometry.jl",
    "title": "PosDefManifold.matP",
    "category": "function",
    "text": "matP(ς::Vector{T}) where T<:RealOrComplex\n\nMatrizize a tangent vector (vector) ς :  vec -> mat.\n\nThis is the function reversing the vecP function,  thus the weighting applied therein is reversed as well.\n\nIf ς=vecP(S) and S is a nn Hermitian matrix,  ς  is a tangent vector of size n(n+1)2.  The result of calling matP(ς) is then nn matrix S.\n\nTo Do: This function needs to be rewritten more efficiently\n\nExamples\n\nusing PosDefManifold\nP=randP(3)\nQ=randP(3)\nG=mean(Fishr, P, Q)\n# projecting P at onto the tangent space at the Fisher mean base point\nS=logMap(Fisher, P, G)\n# vectorize S\nv=vecP(S)\n# Rotate the vector by an orthogonal matrix\nn=Int(size(S, 1)*(size(S, 1)+1)/2)\nU=randP(n)\nz=U*v\n# Get the point in the tangent space\nS=matP(z)\n\n\n\n\n\n"
},

{
    "location": "riemannianGeometry/#Tangent-Space-operations-1",
    "page": "riemannianGeometry.jl",
    "title": "Tangent Space operations",
    "category": "section",
    "text": "Function Description\nlogMap Logarithmic map (from manifold to tangent space)\nexpMap Exponential map (from tangent space to manifold)\nvecP vectorization of matrices in the tangent space\nmatP matrization of matrices in the tangent space (inverse of `vecp)⋅logMap\nexpMap\nvecP\nmatP"
},

{
    "location": "riemannianGeometry/#PosDefManifold.procrustes",
    "page": "riemannianGeometry.jl",
    "title": "PosDefManifold.procrustes",
    "category": "function",
    "text": "procrustes(P::ℍ{T}, Q::ℍ{T}, extremum=\"min\") where T<:RealOrComplex\n\nGiven two positive definite matrices P and Q,  return by default the solution of problem\n\ntextrmargmin_Uδ(PU^*QU),\n\nwhere U varies over the set of unitary matrices and δ() is a  distance or divergence function.\n\nU^*QU is named in physics the unitary orbit of Q.\n\nIf the argument extremum is passed as \"max\", it returns instead the solution of\n\ntextrmargmax_Uδ(PU^*QU).\n\nP and Q must be flagged as Hermitian. See typecasting matrices.\n\nAs it has been shown in Bhatia and Congedo (2019)🎓,  using each of the Fisher, logdet zero, Wasserstein  and the Kullback-Leibler divergence (see logdet α),  the best approximant to P from the unitary orbit of Q  commutes with P and, surprisingly, has the same closed-form expression, namely\n\nU_Q^U_P^* for the argmin and U_Q^U_P^* for the argmax,\n\nwhere U^ denotes the eigenvector matrix of the subscript argument with  eigenvectors in columns sorted by decreasing order of corresponding eigenvalues and  U^ denotes the eigenvector matrix of the subscript argument with  eigenvectors in columns sorted by increasing order of corresponding eigenvalues.\n\nThe same solutions are known since a long time also by solving the extremal  problem here above using the Euclidean metric (Umeyama, 1988).\n\nExamples\n\nusing PosDefManifold\nP=randP(3)\nQ=randP(3)\n# argmin problem\nU=procrustes(P, Q)\n# argmax problem\nV=procrustes(P, Q, \"max\")\n\n\n\n\n\n"
},

{
    "location": "riemannianGeometry/#Procrustes-problems-1",
    "page": "riemannianGeometry.jl",
    "title": "Procrustes problems",
    "category": "section",
    "text": "Function Description\nprocrustes Solution to the Procrustes problem in the manifold of positive definite matrices⋅procrustes"
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
    "text": "This unit contains linear algebra functions useful in relation to the Riemannian  geometry of the manifold of Symmetric Positive Definite (SPD) or  Hermitian Positive Definite (HPD) matrices. In Julia those are Hermitian matrices, see typecasting matrices.In general they take a matrix as input (some may take other arrays as input) and are divided in eight categories depending on what kind of functions thay are and what they give as output:Category Output\n1. Utilities - - -\n2. Matrix normalizations matrix\n3. Boolean functions of matrices matrix\n4. Scalar functions of matrices scalar\n5. Diagonal functions of matrices diagonal matrix\n6. Unitary functions of matrices orthogonal/unitary matrix\n7. Matrix function of matrices matrix\n8. Spectral decompositions of positive matrices spectral function of input\n9. Decompositions involving triangular matrices triangular matrix⋅"
},

{
    "location": "linearAlgebra/#PosDefManifold.typeofMatrix",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.typeofMatrix",
    "category": "function",
    "text": "function typeofMatrix(array::Union{ AnyMatrix,\n                                    AnyMatrixVector,\n                                    AnyMatrixVector₂ })\n\nalias: typeofMat\n\nReturn the type of a matrix, either Hermitian,  Diagonal, LowerTriangular, or Matrix.  Argument array may be a matrix of one of these types, but also one of the following:\n\nℍVector, ℍVector₂, 𝔻Vector, 𝔻Vector₂, 𝕃Vector, 𝕃Vector₂, 𝕄Vector, 𝕄Vector₂.\n\nThose are Array of Matrices types.  See also aliases for the symbols ℍ, 𝔻, 𝕃 and 𝕄.\n\nNote that this function is different from Julia function  typeof,  which returns the concrete type (see example below), thus  cannot be used for typecasting matrices.\n\nExamples\n\nusing LinearAlgebra, PosDefManifold\nP=randP(3) # generate a 3x3 Hermitian matrix\ntypeofMatrix(P) # returns `Hermitian`\ntypeof(P) # returns `Hermitian{Float64,Array{Float64,2}}`\n# typecast P as a `Matrix` M\nM=Matrix(P)\n# typecast M as a matrix of the same type as P and write the result in A\nA=typeofMatrix(P)(M)\n\nPset=randP(3, 4) # generate a set of 4 3x3 Hermitian matrix\n# Pset is an ℍVector type\ntypeofMatrix(Pset) # again returns `Hermitian`\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#PosDefManifold.typeofVector",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.typeofVector",
    "category": "function",
    "text": "function typeofVector(array::Union{ AnyMatrixVector,\n                                    AnyMatrixVector₂ })\n\nalias: typeofVec\n\nReturn the type of a Vector, either HermitianVector,  DiagonalVector, LowerTriangularVector, or MatrixVector.  The aliases of those are, respectvely, ℍVector, 𝔻Vector, 𝕃Vector and 𝕄Vector.  Argument array may be a vector of one of these types, but also one of the following:\n\nℍVector₂, 𝔻Vector₂, 𝕃Vector₂, 𝕄Vector₂.\n\nThose are Array of Matrices types.  See also aliases for the symbols ℍ, 𝔻, 𝕃 and 𝕄.\n\nNote that this function is different from Julia function  typeof  only in that it returns the vector type also for  the ℍVector₂, 𝔻Vector₂, 𝕃Vector₂ and 𝕄Vector₂ types.\n\nExamples\n\nusing LinearAlgebra, PosDefManifold\nP=randP(3, 4) # generate 4 3x3 Hermitian matrix\ntypeofMatrix(P) # returns `Array{Hermitian,1}`\ntypeof(P) # also returns `Array{Hermitian,1}`\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#PosDefManifold.dim",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.dim",
    "category": "function",
    "text": "(1) function dim(X::AnyMatrix, [d])\n(2) function dim(vector::AnyMatrixVector, [d])\n(3) function dim(vector₂::AnyMatrixVector₂, [d])\n\n(1) X is a real or complex Matrix, Diagonal,  LowerTriangular or Hermitian matrix.  Return a 2-tuple containing the dimensions of X,  which is two times the same dimension for all possible types of X  with the exception of the Matrix type, which can be rectangular.  Optionally you can specify a dimension (1 or 2)  to get just the length of that dimension.\n\n(2) vector is an 𝕄Vector, 𝔻Vector, 𝕃Vector or ℍVector type  (see AnyMatrixVector type).  Return a 3-tuple containing the number of matrices it holds  (dimension 1) and their dimensions (dimension 2 and 3).  Optionally you can specify a dimension (1, 2, or 3)  to get just the length of that dimension.\n\n(3) vector₂ is an 𝕄Vector₂, 𝔻Vector₂, 𝕃Vector₂ or ℍVector₂ type  (see AnyMatrixVector type).  Return a 4-tuple containing\n\nthe number of vectors of matrices it holds (dimension 1),\na vector holding the number of matrices in each vector of matrices (dimensions 2),\nthe two dimensions of the matrices (dimension 3 and 4).\n\nOptionally you can specify a dimension (1, 2, 3 or 4)  to get just the length of that dimension.\n\nvector and vector₂ are Array of Matrices types.  See also aliases for the symbols ℍ, 𝔻, 𝕃 and 𝕄.\n\nnote: Nota Bene\nIf you specify a dimension and this is out of the valid range, the function returns zero.Both the vector(2) and the vector₂(3) object are meant to hold matrices living in the same manifold, therefore it is assumed that all matrices they holds are of the same dimension. The dimensions of the matrices are retrived fromthe first matrix in vector(2),\nthe first matrix in the first vector of vector₂(3).\n\nThis function replaces Julia size  function, which cannot be used to retrive dimension for matrix vectors.  It is not possible to overload the size function for matrix vectors  since this causes problems to other Julia functions.\n\nExamples\n\nusing LinearAlgebra, PosDefManifold\n# (1)\nM=randn(3, 4) # generate a 3x4 `Matrix`\ndim(M) # returns (3, 4)\ndim(M, 1) # returns 3\ndim(M, 2) # returns 4\ndim(M, 3) # out of range: returns 0\n\n# (2)\nPset=randP(3, 4) # generate an ℍVector holding 4 3x3 Hermitian matrices\ndim(Pset) # returns (4, 3, 3)\ndim(Pset, 1) # returns 4\ndim(Pset, 2) # returns 3\ndim(Pset, 3) # returns 3\n\n# (3)\n# Generate a set of 4 random 3x3 SPD matrices\nPset=randP(3, 4)\n# Generate a set of 40 random 4x4 SPD matrices\nQset=randP(3, 40)\nA=ℍVector₂([Pset, Qset])\ndim(A) # return (2, [4, 40], 3, 3)\ndim(A, 1) # return 2\ndim(A, 2) # return [4, 40]\ndim(A, 2)[1] # return 4\ndim(A, 3) # return 3\ndim(A, 4) # return 3\ndim(A, 5) # out of range: return 0\n\n# note: to create an ℍVector₂ object holding k ℍVector objects use\nsets=ℍVector₂(undef, k) # and then fill them\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#Utilities-1",
    "page": "linearAlgebra.jl",
    "title": "Utilities",
    "category": "section",
    "text": "Function Description\ntypeofMatrix, typeofMat Return the type of the matrix argument\ntypeofVector, typeofVec Return the type of the matrix vector argument\ndim length of the dimansions of matrices and vectors of matrices⋅typeofMatrix\ntypeofVector\ndim"
},

{
    "location": "linearAlgebra/#PosDefManifold.det1",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.det1",
    "category": "function",
    "text": "function det1(X::AnyMatrix; <tol::Real=0>)\n\nReturn the argument matrix X normalized so as to have unit determinant.  For square positive definite matrices this is the best approximant  from the set of matrices in the special linear group -  see Bhatia and Jain (2014)[🎓].\n\nX can be a real or complex Diagonal, LowerTriangular,  Matrix, or Hermitian matrix. (see AnyMatrix type)\n\nIf the determinant is not greater to tol (which defalts to zero)  a warning is printed and X is returned.\n\nnote: Nota Bene\nThis function is meant for positive definite matrices. Julia may throws an error while computing the determinant if the matrix is defective.\n\nSee Julia det function.\n\nSee also: tr1.\n\nExamples\n\nusing LinearAlgebra, PosDefManifold\nP=randP(5) # generate a random real positive definite matrix 5x5\nQ=det1(P)\ndet(Q) # must be 1\n# using a tolerance\nQ=det1(P; tol=1e-12)\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#PosDefManifold.tr1",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.tr1",
    "category": "function",
    "text": "tr1(X::AnyMatrix)\n\nReturn the argument matrix X normalized so as to have unit trace.\n\nX can be a real or complex Diagonal, LowerTriangular,  Matrix or Hermitian matrix (see AnyMatrix type).\n\nIf the trace is not greater to tol  (which defalts to zero) a warning is printed and X is returned.\n\nSee: Julia trace function.\n\nSee also: tr, det1.\n\nExamples\n\nusing LinearAlgebra, PosDefManifold\nP=randP(5) # generate a random real positive definite matrix 5x5\nQ=tr1(P)\ntr(Q)  # must be 1\n# using a tolerance\nQ=tr1(P; tol=1e-12)\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#PosDefManifold.normalizeCol!",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.normalizeCol!",
    "category": "function",
    "text": "(1) normalizeCol!(X::𝕄{T}, j::Int)\n(2) normalizeCol!(X::𝕄{T}, j::Int, by::Number)\n(3) normalizeCol!(X::𝕄{T}, range::UnitRange)\n(4) normalizeCol!(X::𝕄{T}, range::UnitRange, by::Number)\n                 for all above: where T<:RealOrComplex\n\nGiven a Matrix type X comprised of real or complex elements,\n\n(1) normalize the j^th column to unit norm\n(2) divide the elements of the j^th column by number by\n(3) normalize the columns in range to unit norm\n(4) divide the elements of columns in range  by number by.\n\nby is a number of abstract supertype Number.  It should be an integer, real or complex number.  For efficiency, it should be of the same type as the elements of X.\n\nrange is a UnitRange type.\n\nMethods (1) and (3) call the  BLAS.nrm2  routine for computing the norm of concerned columns.  See Threads.\n\nnote: Nota Bene\nJulia does not allow normalizing the columns of Hermitian matrices. If you want to call this function for an Hermitian matrix see typecasting matrices.\n\nSee norm and also randn  for the example below.\n\nSee also: colNorm, colProd.\n\nExamples\n\nusing PosDefManifold\nX=randn(10, 20)\nnormalizeCol!(X, 2)                  # (1) normalize columns 2\nnormalizeCol!(X, 2, 10.0)            # (2) divide columns 2 by 10.0\nnormalizeCol!(X, 2:4)                # (3) normalize columns 2 to 4\nX=randn(ComplexF64, 10, 20)\nnormalizeCol!(X, 3)                  # (1) normalize columns 3\nnormalizeCol!(X, 3:6, (2.0 + 0.5im)) # (4) divide columns 3 to 5 by (2.0 + 0.5im)\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#Matrix-normalizations-1",
    "page": "linearAlgebra.jl",
    "title": "Matrix normalizations",
    "category": "section",
    "text": "Function Description\ndet1 Normalize the determinant\ntr1 Normalize the trace\nnormalizeCol! Normalize one or more columns⋅det1\ntr1\nnormalizeCol!"
},

{
    "location": "linearAlgebra/#PosDefManifold.ispos",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.ispos",
    "category": "function",
    "text": "    (1) ispos(   λ::Vector{T};\n                <tol::Real=0, rev=true, 🔔=true, msg=\"\">) where T<:Real\n\n    (2) ispos(   Λ::𝔻{T};\n                <tol::Real=0, rev=true, 🔔=true, msg=\"\">) where T<:Real\n\nReturn true if all numbers in (1) real vector λ or in (2) real Diagonal  matrix Λ are not inferior to tol, otherwise return false. This is used,  for example, in spectral functions to check that all eigenvalues are positive.\n\nnote: Nota Bene\ntol defaults to the square root of Base.eps of the type of λ (1)  or Λ (2). This corresponds to requiring positivity beyond about half of  the significant digits.\n\nThe following are <optional keyword arguments>:\n\nIf rev=true the (1) elements in λ or (2) the diagonal elements in Λ will be chacked in reverse order.\n\nThis is done for allowing a very fast check when the elements  are sorted and it is known from where is best to start checking.\n\nIf the result is false:\n\nif =true a bell character will be printed. In most systems this will ring a bell on the computer.\nif string msg is provided, a warning will print msg followed by:\n\n\"at position pos\", where pos is the position where the  first non-positive element has been found.\n\n ## Examples\n using PosDefManifold\n a=[1, 0, 2, 8]\n ispos(a, msg=\"non-positive element found\")\n\n # it will print:\n # ┌ Warning: non-positive element found at position 2\n # └ @ [here julie will point to the line of code issuing the warning]\n\n\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#Boolean-functions-of-matrices-1",
    "page": "linearAlgebra.jl",
    "title": "Boolean functions of matrices",
    "category": "section",
    "text": "Function Description\nispos Check whether a real vector or diagonal matrix are comprised of all positive elementsispos"
},

{
    "location": "linearAlgebra/#PosDefManifold.colProd",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.colProd",
    "category": "function",
    "text": "(1) colProd(X::Union{𝕄{T}, ℍ{T}}, j::Int, l::Int)\n(2) colProd(X::Union{𝕄{T}, ℍ{T}}, Y::Union{𝕄{T}, ℍ{T}}, j::Int, l::Int)\n           for all above: where T<:RealOrComplex\n\n(1) Given a real or complex Matrix or Hermitian matrix X,  return the dot product of the j^th and l^th columns, defined as,\n\nsum_i=1^r big(x_ij^*x_ilbig)\n\nwhere r is the number of rows of X and ^* denotes complex  conjugate (nothing if the matrix is real).\n\n(2) Given real or complex Matrix or Hermitian matrices X and Y,  return the dot product of the j^th column of X and the l^th column  of Y, defined as,\n\nsum_i=1^r big(x_ij^*y_ilbig)\n\nwhere r is the number of rows of X and of Y and ^* is as above.\n\nnote: Nota Bene\nX and of Y may have a different number of columns, but must have the same number of rows.\n\nArguments j and l must be positive integers in range\n\n(1) j,l in 1:size(X, 2),\n(2) j in 1:size(X, 2), l in 1:size(Y, 2).\n\nSee also: normalizeCol!, colNorm.\n\nExamples\n\nusing PosDefManifold\nX=randn(10, 20)\np=colProd(X, 1, 3)\nY=randn(10, 30)\nq=colProd(X, Y, 2, 25)\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#PosDefManifold.sumOfSqr",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.sumOfSqr",
    "category": "function",
    "text": "(1) sumOfSqr(A::Array)\n(2) sumOfSqr(H::ℍ{T})\n(3) sumOfSqr(L::𝕃{T})\n(4) sumOfSqr(D::𝔻{T})\n(5) sumOfSqr(X::Union{𝕄{T}, ℍ{T}}, j::Int)\n(6) sumOfSqr(X::Union{𝕄{T}, ℍ{T}}, range::UnitRange)\n             for (1)-(6) above: where T<:RealOrComplex\n\nalias: ss\n\nReturn\n\n(1) the sum of squares of the elements in an array A of any dimensions.\n(2) as in (1), but for an Hermitian matrix H, using only the lower triangular part.\n(3) as in (1), but for a LowerTriangular matrix L.\n(4) as in (1), but for a Diagonal matrix D (sum of squares of diagonal elements).\n(5) the sum of square of the j^th column of a Matrix or Hermitian X.\n(6) the sum of square of the columns of a Matrix or Hermitian X in a given range.\n\nAll methods support real and complex matrices.\n\nOnly method (1) works for arrays of any dimensions.\n\nMethods (1)-(4) return the square of the Frobenius norm.\n\nFor method (5), j is a positive integer in range 1:size(X, 1).\n\nFor method (6), range is a UnitRange type.\n\nSee also: colNorm, sumOfSqrDiag, sumOfSqrTril.\n\nExamples\n\nusing PosDefManifold\nX=randn(10, 20)\nsum2=sumOfSqr(X)        # (1) sum of squares of all elements\nsum2=sumOfSqr(X, 1)     # (2) sum of squares of elements in column 1\nsum2=sumOfSqr(X, 2:4)   # (3) sum of squares of elements in column 2 to 4\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#PosDefManifold.sumOfSqrDiag",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.sumOfSqrDiag",
    "category": "function",
    "text": "sumOfSqrDiag(X::AnyMatrix)\n\nalias: ssd\n\nSum of squares of the diagonal elements in real or complex Matrix,  Diagonal, Hermitian or LowerTriangular matrix X.  If X is rectangular (which can be only if it is of the Matrix type),  the main diagonal is considered.\n\nSee AnyMatrix type\n\nSee also: sumOfSqr, sumOfSqrTril.\n\nExamples\n\nusing LinearAlgebra, PosDefManifold\nX=randn(10, 20)\nsumDiag2=sumOfSqrDiag(X) # (1)\nsumDiag2=sumOfSqrDiag(𝔻(X)) # (2) 𝔻=LinearAlgebra.Diagonal\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#PosDefManifold.colNorm",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.colNorm",
    "category": "function",
    "text": "colNorm(X::Union{𝕄{T}, ℍ{T}}, j::Int) where T<:RealOrComplex\n\nGiven a real or complex Matrix or Hermitian matrix X,  return the Euclidean norm of its j^th column.\n\nThis function calls the  BLAS.nrm2  routine. See Threads.\n\nSee also: normalizeCol!, colProd, sumOfSqr.\n\nExamples\n\nusing PosDefManifold\nX=randn(10, 20)\nnormOfSecondColumn=colNorm(X, 2)\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#PosDefManifold.sumOfSqrTril",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.sumOfSqrTril",
    "category": "function",
    "text": "sumOfSqrTril(X::AnyMatrix, k::Int=0)\n\nalias: sst\n\nGiven a real or complex Matrix, Diagonal, Hermitian or  LowerTriangular matrix X (see AnyMatrix type),  return the sum of squares of the elements  in its lower triangle up to the k^th underdiagonal.\n\nMatrix X may be rectangular.\n\nk must be in range\n\n1-size(X, 1):c-1 for X Matrix, Diagonal or Hermitian,\n1-size(X, 1):0 for X LowerTriangular.\n\nFor X Diagonal the result is\n\n0 if k0,\nthe sum of the squares of the diagonal elements otherwise.\n\nSee julia tril(M, k::Integer) function  for numbering of diagonals.\n\nSee also: sumOfSqr, sumOfSqrDiag.\n\nExamples\n\nusing PosDefManifold\nA=[4. 3.; 2. 5.; 1. 2.]\n#3×2 Array{Float64,2}:\n# 4.0  3.0\n# 2.0  5.0\n# 1.0  2.0\n\ns=sumOfSqrTril(A, -1)\n# 9.0 = 1²+2²+2²\n\ns=sumOfSqrTril(A, 0)\n# 50.0 = 1²+2²+2²+4²+5²\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#LinearAlgebra.tr",
    "page": "linearAlgebra.jl",
    "title": "LinearAlgebra.tr",
    "category": "function",
    "text": "(1) tr(P::ℍ{T}, Q::ℍ{T})\n(2) tr(P::ℍ{T}, M::𝕄{T})\n(3) tr(D::𝔻{T}, H::Union{ℍ{T}, 𝕄{T}})\n(4) tr(H::Union{ℍ{T}, 𝕄{T}}, D::𝔻{T})\n      for all above: where T<:RealOrComplex\n\nGiven (1) two Hermitian positive definite matrix P and Q,  return the trace of the product PQ.  This is real even if P and Q are complex.\n\nP must always be flagged as Hermitian. See typecasting matrices.\n\nIn (2) Q is a Matrix object,  in which case return\n\na real trace if the product PQ is real or if it has all positive real eigenvalues.\na complex trace if the product PQ is not real and has complex eigenvalues.\n\nMethods (3) and (4) return the trace of the product DH or HD,  where D is a Diagonal matrix and H an Hermitian  or Matrix object. The result is of the same type as the input matrices.\n\nFor all methods all arguments must be of the same type.\n\nMath\n\nLet P and Q be Hermitian matrices, using the properties of the trace  (e.g., the cyclic property and the similarity invariance) you can use this  function to fast compute the trace of several expressions. For example:\n\ntextrmtr(PQ)=textrmtr(P^12QP^12)\n\nand\n\ntextrmtr(PQP)=textrmtr(P^2Q) (see example below).\n\nSee: trace.\n\nSee also: DiagOfProd, tr1.\n\nExamples\n\nusing PosDefManifold\nP=randP(ComplexF64, 5) # generate a random complex positive definite matrix 5x5\nQ=randP(ComplexF64, 5) # generate a random complex positive definite matrix 5x5\ntr(P, Q) ≈ tr(P*Q) ? println(\" ⭐ \") : println(\" ⛔ \")\ntr(P, Q) ≈ tr(sqrt(P)*Q*sqrt(P)) ? println(\" ⭐ \") : println(\" ⛔ \")\ntr(sqr(P), Q) ≈ tr(P*Q*P) ? println(\" ⭐ \") : println(\" ⛔ \")\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#PosDefManifold.quadraticForm",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.quadraticForm",
    "category": "function",
    "text": "(1) quadraticForm(v::Vector{T}, P::ℍ{T}) where T<:Real\n(2) quadraticForm(v::Vector{T}, L::𝕃{T}) where T<:Real\n(3) quadraticForm(v::Vector{T}, X::𝕄{T}, forceLower::Bool=false) where T<:Real\n(4) quadraticForm(v::Vector{S}, X::Union{𝕄{S}, ℍ{S}}) where S<:Complex\n\nalias: qf\n\n(1) Given a real vector v and a real Hermitian matrix P,  compute the quadratic form\n\nv^TPv,\n\nwhere the superscript T denotes transpose,  using only the lower triangular part of P.\n\n(2) As in (1), given a real vector v  and the LowerTriangular view L of a real symmetric matrix.\n\n(3) As in (1), given a real vector v  and a real generic Matrix M, if forceLower=true. If forceLower=false,  the product v^TMv is evaluated using the whole matrix M.\n\n(4) Quadratic form v^HPv, where superscript H denotes complex conjugate  and transpose, for a complex vector v and  Matrix or Hermitian matrix.  The whole matrix is then used.\n\nMath\n\nFor v and X real and X symmetric, the quadratic form is\n\nsum_i(v_i^2x_ii)+sum_ij(2v_iv_jx_ij).\n\nThis is used in (1), (2) and (3).\n\nExamples\n\nusing PosDefManifold\nP=randP(5) # generate a random real positive definite matrix 5x5\nv=randn(5)\nq1=quadraticForm(v, P) # or q1=quad(v, P)\n# obtain a lower Triangular view of P\nL=LowerTriangular(Matrix(P)) # or L=𝕃(Matrix(P))\nq2=quadraticForm(v, L)\nq1 ≈ q2 ? println(\" ⭐ \") : println(\" ⛔ \")\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#PosDefManifold.fidelity",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.fidelity",
    "category": "function",
    "text": "fidelity(P::ℍ{T}, Q::ℍ{T}) where T<:RealOrComplex\n\nGiven two positive definte Hermitian matrices P and Q,  return their fidelity:\n\ntrbig(P^12QP^12big)^12\n\nThis is used in quantum physics and is related to the  Wasserstein metric. See for example Bhatia, Jain and Lim (2019b)🎓.\n\nExamples\n\nusing PosDefManifold\nP=randP(5);\nQ=randP(5);\nf=fidelity(P, Q)\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#Scalar-functions-of-matrices-1",
    "page": "linearAlgebra.jl",
    "title": "Scalar functions of matrices",
    "category": "section",
    "text": "Function Description\ncolProd Sum of products of the elements in two columns\nsumOfSqr, ss Sum of squares of all elements or of specified columns\nsumOfSqrDiag, ssd Sum of squares of the diagonal elements\ncolNorm Eucliden norm of a column\nsumOfSqrTril, sst Sum of squares of the lower triangle elements up to a given underdiagonal\ntr Fast trace of the product of two Hermitian matrices\nquadraticForm, qf Fast quadratic form\nfidelity (Quantum) Fidelity of two positive matrices⋅colProd\nsumOfSqr\nsumOfSqrDiag\ncolNorm\nsumOfSqrTril\ntr\nquadraticForm\nfidelity"
},

{
    "location": "linearAlgebra/#PosDefManifold.fDiag",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.fDiag",
    "category": "function",
    "text": "fDiag(func::Function, X::AnyMatrix, k::Int=0)\n\nalias: 𝑓𝔻\n\nApplies function func element-wise to the elements of the k^th  diagonal of real or complex Diagonal, LowerTriangular, Matrix  or Hermitian matrix X and return a diagonal matrix with these elements.  X must be square in all cases, but for the 𝕄=Matrix type argument,  in which case it may be of dimension r⋅c, with r ≠ c.\n\nSee julia tril(M, k::Integer) function  for numbering of diagonals.\n\nBt default the main diagonal is considered.\n\nIf X is Diagonal, k must be zero (main diagonal).\nIf X is LowerTriangular, k cannot be positive.\n\nNote that if X is rectangular the dimension of the result depends  on the size of X and on the chosen diagonal.  For example,\n\nr ≠ c and k=0 (main diagonal), the result will be of dimension min(r,c)⋅min(r,c),\nX 3⋅4 and k=-1, the result will be 2⋅2,\nX 3⋅4 and k=1, the result will be 3⋅3, etc.\n\nnote: Nota Bene\nThe function func must support the func. syntax and therefore must be able to apply element-wise to the elements of the chosen diagonal (this includes anonymous functions).\n\nIf the input matrix is complex, the function `func`\nmust be able to support complex arguments.\n\nSee also: DiagOfProd, tr.\n\nExamples\n\nusing PosDefManifold\nP=randP(5) # use P=randP(ComplexF64, 5) for generating an Hermitian matrix\nD=fDiag(inv, P, -1)   # diagonal matrix with the inverse of the first sub-diagonal of P\n(Λ, U) = evd(P)       # Λ holds the eigenvalues of P, see evd\nΔ=fDiag(log, Λ)       # diagonal matrix with the log of the eigenvalues\nΔ=fDiag(x->x^2, Λ)    # using an anonymous function for the square of the eigenvalues\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#PosDefManifold.DiagOfProd",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.DiagOfProd",
    "category": "function",
    "text": "DiagOfProd(P::ℍ{T}, Q::ℍ{T}) where T<:RealOrComplex\n\nalias: dop\n\nReturn the Diagonal matrix holding the diagonal of the product PQ  of two Hermitian matrices P and Q. Only the diagoanl part  of the product is computed.\n\nSee also: tr, fDiag.\n\nExamples\n\nusing PosDefManifold, LinearAlgebra\nP, Q=randP(5), randP(5)\nDiagOfProd(P, Q)≈Diagonal(P*Q) ? println(\"⭐ \") : println(\"⛔ \")\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#Diagonal-functions-of-matrices-1",
    "page": "linearAlgebra.jl",
    "title": "Diagonal functions of matrices",
    "category": "section",
    "text": "Function Description\nfDiag, 𝑓𝔻 Elemen-wise functions of matrix diagonals\nDiagOfProd, dop Diagonal of the product of two matrices⋅fDiag\nDiagOfProd"
},

{
    "location": "linearAlgebra/#PosDefManifold.mgs",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.mgs",
    "category": "function",
    "text": "mgs(X::𝕄{T}, numCol::Int=0) where T<:RealOrComplex\n\nModified (stabilized) Gram-Schmidt orthogonalization  of the columns of square or tall matrix X, which can be comprised of real  or complex elements.  The orthogonalized X is returned by the function. X is not changed.\n\nAll columns are orthogonalized by default. If instead argument numCol is provided,  then only the first numCol columns of X are orthogonalized.  In this case only the firt numCol columns will be returned.\n\nExamples\n\nusing LinearAlgebra, PosDefManifold\nX=randn(10, 10);\nU=mgs(X)        # result is 10⋅10\nU=mgs(X, 3)     # result is 10⋅3\nU\'*U ≈ I ? println(\" ⭐ \") : println(\" ⛔ \")\n# julia undertands also:\nU\'U ≈ I ? println(\" ⭐ \") : println(\" ⛔ \")\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#Unitary-functions-of-matrices-1",
    "page": "linearAlgebra.jl",
    "title": "Unitary functions of matrices",
    "category": "section",
    "text": "Function Description\nmgs Modified Gram-Schmidt orthogonalization⋅mgs"
},

{
    "location": "linearAlgebra/#PosDefManifold.fVec",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.fVec",
    "category": "function",
    "text": "(1) fVec(f::Function, 𝐏::AnyMatrixVector;\n	 	<w::Vector=[], ✓w=false, allocs=[]>)\n\n(2) fVec(f::Function, g::Function, 𝐏::AnyMatrixVector;\n		<w::Vector=[], ✓w=false, allocs=[]>)\n\nGiven a 1d array 𝐏=P_1P_k of k matrices  of the 𝕄Vector type, 𝔻Vector type, 𝕃Vector type or  ℍVector type and an optional non-negative real weights vector  w=w_1w_k, return expression\n\n(1)hspace6ptf_i=1^k(w_iP_i),\n\nor\n\n(2)hspace6ptf_i=1^k(w_ig(P_i)),\n\nwhere f is a linear matrix function iterating over all elements of 𝐏,  typically the mean or sum function and g is whatever matrix function  applying to each matrix P_k, such as exp, log,sqrt`, etc,  and anonymous functions.\n\nThis function is always multi-threaded. It works by partitioning the k  operations required by the f function in several groups,  passing each group to a separate thread and applying the f function  again on the intermediate results (that\'s why f must be linear).  This function allows a gain in computational time only when the number of  matrices (1) and/or their size (2) is high. Use mean and sum otherwise.  For the number of threads Julia is instructed to use see Threads.\n\n<optional keword argument> allocs allows to pass pre-allocated memory  for holding the intermediate result of each thread.  Argument allocs must be a vector of as many matrices as threads and where  the matrices have the same dimension as the the matrices in 𝐏  (see the example here below). Using this option is worthwhile only  if the size of the matrices is very high and/or when fVec is to be  called repeatedly on many vector of matrices, where the matrices  have always the same size, so that one allocation works for all calls.\n\nIf <optional keyword argument> ✓w=true is passed, the weights are  normalized so as to sum up to 1, otherwise they are used as they are passed.  This option is provided to allow calling this function repeatedly without  normalizing the same weights vector each time. By default ✓w is false.\n\nnote: Nota Bene\nContrarily to Julia mean and sum function (v 1.1.0) the fVec function returns a matrix of the same type of the matrices in 𝐏.\n\nExamples\n\nusing LinearAlgebra, PosDefManifold\nPset=randP(4, 1000); # generate 1000 positive definite 4x4 matrices\nmean(Pset) # arithmetic mean calling Julia function\nThreads.nthreads() # check that at least two threads are available\nfVec(mean, Pset) # multi-threaded arithmetic mean\n\ninv(mean(inv, Pset)) # Harmonic mean calling Julia function\ninv(fVec(mean, inv, Pset)) # multi-threaded Harmonic mean\n\nexp(mean(log, Pset)) # log Euclidean mean calling Julia function\nexp(fVec(mean, log, Pset)) # multi-threaded log Euclidean mean\n\n# notice that Julia `exp` function has changed the type of the result\n# to `Symmetric`. To obtain an `Hermitian` output use\nℍ(exp(fVec(mean, log, Pset)))\n\nw=(randn(1000)).^2\nw=w./sum(w)  		# generate normalized random weights\n\n# weighted arithmetic mean calling Julia function\nsum(Pset[i]*w[i] for i=1:length(w))\n# multi-threaded weighted arithmetic mean\nfVec(sum, Pset, w=w)\n\n# weighted harmonic mean calling Julia function\ninv(sum(inv(Pset[i])*w[i] for i=1:length(w)))\n# multi-threaded weighted harmonic mean\ninv(fVec(sum, inv, Pset, w=w))\n\n# pre-allocating memory\nPset=randP(100, 1000); # generate 1000 positive definite 100x100 matrices\nQset=MatrixVector(repeat([Pset[1]], nthreads()))\nfVec(mean, log, Pset, allocs=Qset)\n\n# How much computing time we save ?\n# (example min time obtained with 4 threads & 4 BLAS threads)\nusing BenchmarkTools\n# standard Julia function\n@benchmark(mean(log, Pset)) 					# (5.271 s)\n# fVec\n@benchmark(fVec(mean, log, Pset))				# (1.540 s)\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#PosDefManifold.congruence",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.congruence",
    "category": "function",
    "text": "(1) congruence(B::AnyMatrix, P::ℍ)\n(2) congruence(B::AnyMatrix, 𝐏::ℍVector)\n\nalias: cong\n\n(1) Return as an Hermitian matrix the congruent transformation  of Hermitian matrix P\n\nBPB^H,\n\nfor B Hermitian, LowerTriangular, Diagonal or a general Matrix.\n\n(2) Return an ℍVector type holding the congruent transformations\n\nBP_kB^H,\n\nfor all k matrices in ℍVector type 𝐏=P_1P_k.\n\nMethod (2) is multi-threaded. See Threads.\n\nExamples\n\nusing LinearAlgebra, PosDefManifold\n\n# (1)\nP=randP(3) # generate a 3x3 positive matrix\nM=randn(3, 3)\nC=congruence(M, P) # = M*P*M\'\n\n# (2)\nPset=randP(4, 100); # generate 100 positive definite 4x4 matrices\nM=randn(4, 4)\nQset=cong(M, Pset) # = [M*Pset_1*M\',...,M*Pset_k*M\'] as an ℍVector type\n\n# recenter the matrices in Pset to their Fisher mean:\nQset=cong(invsqrt(mean(Fisher, Pset; ⏩=true)), Pset)\n\n# as a check, the Fisher mean of Qset is now the identity\nmean(Fisher, Qset; ⏩=true)≈I ? println(\"⭐\") : println(\"⛔\")\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#Matrix-function-of-matrices-1",
    "page": "linearAlgebra.jl",
    "title": "Matrix function of matrices",
    "category": "section",
    "text": "Function Description\nfVec General function for multi-threaded computation of means and sums of matrix vectors\ncongruence, cong Compute congruent transformations⋅fVec\ncongruence"
},

{
    "location": "linearAlgebra/#PosDefManifold.evd",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.evd",
    "category": "function",
    "text": "evd(S::ℍ{T}) where T<:RealOrComplex\n\nGiven a positive semi-definite matrix S,  returns a 2-tuple (Λ U), where U is the matrix holding in columns  the eigenvectors and Λ is the matrix holding the eigenvalues on the diagonal.  This is the output of Julia  eigen  function in UΛU=S form.\n\nAs for the eigen function, the eigenvalues and associated  eigenvectors are sorted by increasing values of eigenvalues.\n\nS may be real or complex and must be flagged by Julia as Hermitian.  See typecasting matrices.\n\nSee also: spectralFunctions.\n\nExamples\n\nusing PosDefManifold\nA=randn(3, 3);\nS=ℍ(A+A\');\nΛ, U=evd(S); # which is equivalent to (Λ, U)=evd(P)\n(U*Λ*U\') ≈ S ? println(\" ⭐ \") : println(\" ⛔ \")\n# => UΛU\'=S, UΛ=SU, ΛU\'=U\'S\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#PosDefManifold.spectralFunctions",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.spectralFunctions",
    "category": "function",
    "text": "(1) spectralFunctions(P::ℍ{T}, func) where T<:RealOrComplex\n(2) spectralFunctions(D::𝔻{S}, func) where S<:Real\n\n(1) This is the mother function for all spectral functions of eigenvalues implemented  in this library, which are:\n\npow     (power),\nisqrt   (inverse square root).\n\nThe function sqr (square) does not use it, as it can be obtained more  efficiently by simple multiplication.\n\nYou can use this function if you need another spectral function of eigenvalues  besides those and those already implemented in the standard package LinearAlgebra.  In general, you won\'t call it directly.\n\nfunc is the function that will be applied on the eigenvalues.\n\nP must be flagged as Hermitian. See typecasting matrices.  It must be a positive definite or positive semi-definite matrix,  depending on func.\n\nA special method is provided for real Diagonal matrices (2).\n\nnote: Nota Bene\nThe function func must support the func. syntax and therefore must be able to apply element-wise to the eigenvalues (those include anonymous functions).\n\nMaths\n\nThe definition of spectral functions for a positive definite matrix P  is at it follows:\n\nfbig(Pbig)=Ufbig(Λbig)U^H\n\nwhere U is the matrix holding in columns the eigenvectors of P,  Λ is the matrix holding on diagonal its eigenvalues and f is  a function applying element-wise to the eigenvalues.\n\nSee also: evd.\n\nExamples\n\nusing LinearAlgebra, PosDefManifold\nn=5\nP=randP(n) # P=randP(ComplexF64, 5) to generate an Hermitian complex matrix\nnoise=0.1;\nQ=spectralFunctions(P, x->x+noise) # add white noise to the eigenvalues\ntr(Q)-tr(P) ≈ noise*n ? println(\" ⭐ \") : println(\" ⛔ \")\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#PosDefManifold.pow",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.pow",
    "category": "function",
    "text": "(1) pow(P::ℍ{T}, args...) where T<:RealOrComplex\n(2) pow(D::𝔻{S}, args...) where S<:Real\n\n(1) Given a positive semi-definite Hermitian matrix P, return the power  P^r_1 P^r_2  for any number of exponents r_1 r_2.  It returns a tuple comprising as many elements as arguments passed after P.\n\nP must be flagged as Hermitian. See typecasting matrices.\n\narg1 arg2 are real numbers.\n\nA special method is provided for real Diagonal matrices (2).\n\nSee also: invsqrt.\n\nExamples\n\nusing LinearAlgebra, PosDefManifold\nP=randP(5);     # use P=randP(ComplexF64, 5) for generating an Hermitian matrix\nQ=pow(P, 0.5);            # =>  QQ=P\nQ, W=pow(P, 0.5, -0.5);\nW*P*W ≈ I ? println(\" ⭐ \") : println(\" ⛔ \")\nQ*Q ≈ P ? println(\" ⭐ \") : println(\" ⛔ \")\nR, S=pow(P, 0.3, 0.7);\nR*S ≈ P ? println(\" ⭐ \") : println(\" ⛔ \")\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#PosDefManifold.invsqrt",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.invsqrt",
    "category": "function",
    "text": "(1) invsqrt(P{T}::ℍ) where T<:RealOrComplex\n(2) invsqrt(D{S}::𝔻) where S<:Real\n\nGiven a positive definite Hermitian matrix P,  compute the inverse of the principal  square root P^-12.\n\nP must be flagged as Hermitian. See typecasting matrices.\n\nA special method is provided for real Diagonal matrices (2).\n\nMaths\n\nThe principal square root of a positive definite matrix P is the only  symmetric (if P is real) or Hermitian (if P is complex) square root.  Its inverse P^-12 is also named the whitening or sphering  matrix sinceP^-12PP^-12=I.\n\nSee: typecasting matrices.\n\nSee also: pow.\n\nExamples\n\nusing LinearAlgebra, PosDefManifold\nP=randP(ComplexF64, 5);\nQ=invsqrt(P);\nQ*P*Q ≈ I ? println(\" ⭐ \") : println(\" ⛔ \")\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#PosDefManifold.sqr",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.sqr",
    "category": "function",
    "text": "(1) sqr(P::ℍ{T}) where T<:RealOrComplex\n(2) sqr(X::Union{𝕄{T}, 𝕃{T}, 𝔻{S}}) where T<:RealOrComplex where S<:Real\n\n(1) Given a positive semi-definite Hermitian matrix P,  compute its square P^2.\n\nP must be flagged as Hermitian. See typecasting matrices.\n\nA method is provided also for generic matrices of the Matrix type,  LowerTriangular matrices and real Diagonal matrices (2). The output  is of the same type as the input.\n\nSee also: pow.\n\nExamples\n\nusing PosDefManifold\nP=randP(5);\nP²=sqr(P);  # =>  P²=PP\nsqrt(P²)≈ P ? println(\" ⭐ \") : println(\" ⛔ \")\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#PosDefManifold.powerIterations",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.powerIterations",
    "category": "function",
    "text": "powerIterations(H::Union{ℍ{T}, 𝕄{T}}, q::Int;\n<evalues=false, tol::Real=0, maxiter=300, ⍰=false>)  where T<:RealOrComplex\n\npowerIterations(L::𝕃{S}, q::Int;\n<evalues=false, tol::Real=0, maxiter=300, ⍰=false)> where S<:Real\n\nalias: powIter\n\n(1) Compute the q eigenvectors associated to the q largest (real) eigenvalues  of real or complex Hermitian or Matrix H using the  power iterations +  Gram-Schmidt orthogonalization as suggested by Strang.  The eigenvectors are returned with the same type as the elements of H.\n\nH must have real eigenvalues, that is, it must be a symmetric matrix if it is real  or an Hermitian matrix if it is complex.\n\n(2) as in (1), but using only the LowerTriangular view L of a matrix.  This option is available only for real matrices (see below).\n\nThe following are <optional keyword arguments>:\n\ntol is the tolerance for the convergence of the power method (see below),\nmaxiter is the maximum number of iterations allowed for the power method,\nif =true, the convergence of all iterations will be printed,\nif evalues=true, return the 4-tuple (Λ U iterations covergence),\nif evalues=false return the 3-tuple (U iterations covergence).\n\nnote: Nota Bene\nDifferently from the evd function, the eigenvectors and eigenvalues are sorted by decreasing order of eigenvalues.If H is Hermitian and real, only its lower triangular part is used for computing the power iterations, like in (2). In this case the BLAS.symm routine is used. Otherwise the BLAS.gemm routine is used. See Threads.tol defaults to 100 times the square root of Base.eps of the type of H. This corresponds to requiring the relative convergence criterion over two successive iterations to vanish for about half the significant digits minus 2.\n\nSee also: mgs.\n\nExamples\n\nusing LinearAlgebra, PosDefManifold\n# Generate an Hermitian (complex) matrix\nH=randP(ComplexF64, 10);\n# 3 eigenvectors and eigenvalues\nU, iterations, convergence=powIter(H, 3, ⍰=true)\n# all eigenvectors\nΛ, U, iterations, convergence=powIter(H, size(H, 2), evalues=true, ⍰=true);\nU\'*U ≈ I && U*Λ*U\'≈H ? println(\" ⭐ \") : println(\" ⛔ \")\n\n# passing a `Matrix` object\nΛ, U, iterations, convergence=powIter(Matrix(H), 3, evalues=true)\n\n# passing a `LowerTriangular` object (must be a real matrix in this case)\nL=𝕃(randP(10))\nΛ, U, iterations, convergence=powIter(L, 3, evalues=true)\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#Spectral-decompositions-of-positive-matrices-1",
    "page": "linearAlgebra.jl",
    "title": "Spectral decompositions of positive matrices",
    "category": "section",
    "text": "Function Description\nevd Eigenvalue-Eigenvector decomposition of a matrix in UΛU=P form\nspectralFunctions Mother function for creating spectral functions of eigenvalues\npow Power of a positive matrix for any number of exponents in one pass\ninvsqrt Principal square root inverse (whitening) of a positive matrix\nsqr Square of a positive matrix\npowerIterations, powIter Power method for estimating any number of eigenvectors and associated eigenvalues⋅evd\nspectralFunctions\npow\ninvsqrt\nsqr\npowerIterations"
},

{
    "location": "linearAlgebra/#PosDefManifold.choL",
    "page": "linearAlgebra.jl",
    "title": "PosDefManifold.choL",
    "category": "function",
    "text": "(1) choL(P::ℍ{T}) where T<:RealOrComplex\n(2) choL(D::𝔻{S}) where S<:Real\n\n(1) Given a real or complex positive definite Hermitian matrix P,  return the Cholesky lower triangular factor L  such that LL^H=P. To obtain L^H or both L and L^H, use instead  julia function cholesky.\n\nOn output, L is of type LowerTriangular.\n\n(2) For a real Diagonal matrix D, return D^12.\n\nExamples\n\nusing PosDefManifold\nP=randP(5);\nL=choL(P);\nL*L\'≈ P ? println(\" ⭐ \") : println(\" ⛔ \")\n\n\n\n\n\n"
},

{
    "location": "linearAlgebra/#Decompositions-involving-triangular-matrices-1",
    "page": "linearAlgebra.jl",
    "title": "Decompositions involving triangular matrices",
    "category": "section",
    "text": "Function Description\nchoL Lower triangula factor of Cholesky decomposition⋅choL"
},

{
    "location": "signalProcessing/#",
    "page": "signalProcessing.jl",
    "title": "signalProcessing.jl",
    "category": "page",
    "text": ""
},

{
    "location": "signalProcessing/#PosDefManifold.randChi²",
    "page": "signalProcessing.jl",
    "title": "PosDefManifold.randChi²",
    "category": "function",
    "text": "randChi²(df::Int)\n\nalias: randχ²\n\nGenerate a random variable distributed as a chi-squared with df  degrees of freedom.\n\nIt uses the Wilson–Hilferty transformation for df>=20 -  see chi-squared distribution.\n\nExamples\n\nusing Plots, PosDefManifold\nchi=[randχ²(2) for i=1:10000]\nhistogram(chi) # needs Plots package. Check your plots back-end.\n\n\n\n\n\n"
},

{
    "location": "signalProcessing/#PosDefManifold.randEigvals",
    "page": "signalProcessing.jl",
    "title": "PosDefManifold.randEigvals",
    "category": "function",
    "text": "randEigvals(n::Int; <df::Int=2, eigvalsSNR::Real=10e3>)\n\nalias: randλ\n\nGenerate an n-vector of random real positive eigenvalues.  The eigenvalues are generated as in function randΛ(randEigvalsMat),  the syntax of which is used.\n\nSee also: randU (randUnitaryMat), randP (randPosDefMat).\n\nExamples\n\nusing Plots, PosDefManifold\nλ=sort(randλ(10), rev=true)\nσ=sort(randλ(10, eigvalsSNR=10), rev=true)\nplot(λ) # needs Plots package. Check your plots back-end.\nplot!(σ) # needs Plots package. Check your plots back-end.\n\n\n\n\n\n"
},

{
    "location": "signalProcessing/#PosDefManifold.randEigvalsMat",
    "page": "signalProcessing.jl",
    "title": "PosDefManifold.randEigvalsMat",
    "category": "function",
    "text": "(1) randEigvalsMat(n::Int; <df::Int=2, eigvalsSNR::Real=10e3>)\n(2) randEigvalsMat(n::Int, k::Int; <df::Int=2, eigvalsSNR::Real=10e3>)\n\nalias: randΛ\n\n(1) Generate an nn diagonal matrix of random real positive eigenvalues.  (2) An array 1d (of 𝔻Vector type) of k matrices of the kind in (1)\n\nThe eigenvalues are generated according to model\n\nλ_i=χ_df^2+ηhspace6pttextrmforhspace2pti=1n\n\nwhere\n\nχ_df^2 (signal term) is randomly distributed as a chi-square with df degrees of freedom,\nη is a white noise term, function of <keyword argument> eigvalsSNR, such that\n\ntextrmeigenvalues SNR=mathbbEbig(sum_i=1^nλ_ibig)bignη\n\nThe expected sum mathbbEbig(sum_i=1^nλ_ibig) here above is the  expected variance of the signal term, i.e., n(df), since the expectation  of a random chi-squared variable is equal to its degrees of freedom.\n\nIf eigvalsSNR=Inf is passed as argument, then η is set to zero, i.e.,  no white noise is added. In any case eigvalsSNR must be positive.\n\nNote that with the default value of <keyword argument> df (df=2)  the generating model assumes that the eigenvalues  have exponentially decaying variance, which is often observed on real data.\n\nnote: Nota Bene\nThe <keyword argument> eigvalsSNR expresses the expected eigenvalues SNR (signal-to-noise ratio), not the real one, and is not expressed in decibels, but as the expected SNR variance ratio.\n\nThis function is used by function randP (randPosDefMat) to generate  random positive definite matrices with added white noise in order  to emulate eigenvalues observed in real data and to  improve the conditioning of the generated matrices with respect to inversion.\n\nSee also: randλ (randEigvals), randU (randUnitaryMat),  randP (randPosDefMat), randχ² (randChi²).\n\nExamples\n\nusing PosDefManifold\n# (1)\nn=3;\nU=randU(n);\nΛ=randΛ(n, eigvalsSNR=100)\nP=U*Λ*U\' # generate an SPD matrix\nusing LinearAlgebra\nQ=ℍ(U*Λ*U\') # generate an SPD matrix and flag it as \'Hermitian\'\n\n# (2) generate an array of 10 matrices of simulated eigenvalues\nDvec=randΛ(n, 10)\n\n\n\n\n\n"
},

{
    "location": "signalProcessing/#PosDefManifold.randUnitaryMat",
    "page": "signalProcessing.jl",
    "title": "PosDefManifold.randUnitaryMat",
    "category": "function",
    "text": "(1) randUnitaryMat(n::Int)\n(2) randUnitaryMat(::Type{Complex{T}}, n::Int)\n\naliases: randOrthMat, randU\n\nGenerate a random nn\n\n(1) orthogonal matrix (real)\n(2) unitary matrix (complex)\n\nThe matrices are generated running the modified (stabilized)  Gram-Schmidt orthogonalization  procedure (mgs) on an nn matrix filled with random Gaussian elements.\n\nSee also: randΛ (randEigvals), randP (randPosDefMat).\n\nExamples\n\nusing PosDefManifold\nn=3;\nX=randU(n)*sqrt(randΛ(n))*randU(n)\'  # (1) generate a random square real matrix\n\nU=randU(ComplexF64, n);\nV=randU(ComplexF64, n);\nY=U*sqrt(randΛ(n))*V\' # (2) generate a random square complex matrix\n\n\n\n\n\n"
},

{
    "location": "signalProcessing/#PosDefManifold.randPosDefMat",
    "page": "signalProcessing.jl",
    "title": "PosDefManifold.randPosDefMat",
    "category": "function",
    "text": "(1) randPosDefMat(n::Int; <df::Int=2, eigvalsSNR::Real=10e3>)\n(2) randPosDefMat(::Type{Complex{T}}, arguments in (1))\n(3) randPosDefMat(n::Int, k::Int; df::Int=2, eigvalsSNR::Real=10e3, SNR::Real=100)\n(4) randPosDefMat(::Type{Complex{T}}, arguments in (3))\n\nalias: randP\n\nGenerate\n\n(1) one random Hermitian positive definite matrix (real) of size nn\n(2) one random Hermitian positive definite matrix (complex) of size nn\n(3) an array 1d (of ℍVector type) of k matrices of the kind in (1)\n(4) an array 1d (of ℍVector type) of k matrices of the kind in (2).\n\nFor (1) and (2) the matrix is generated according to model\n\nUΛU^H+ηI,\n\nwhere U is a random orthogonal (1) or unitary (2) matrix generated by  function randU(randUnitaryMat) and Λ, η are a positive definite  diagonal matrix and a non-negative scalar depending on <keywords arguments>  df and eigvalsSNR randomly generated calling function  randΛ(randEigvalsMat).\n\nFor (3) and (4) the k matrices are generated according to model\n\n(UΛ_iU^H+ηI)+φ(V_iΔ_iV_i^H+ηI)hspace8pt  Eq.[1]\n\nwhere\n\nU and the V_i are random (3) orthogonal/(4) unitary matrices,\nΛ_i and Δ_i are positive definite diagonal matrices\nη is a non-negative scalar.\n\nAll variables here above are randomly generated as in (1) and (2)\n\nφ is adjusted so as to obtain a desired output SNR (keyword argument) (signal-to-noise ratio), such as\n\nSNR=fracdisplaystylesum_i=1^ktextrmtr(UΛ_kU^H+ηI)displaystylesum_i=1^ktextrmtrφ(UΔ_kU^H+ηI).\n\nnote: Nota Bene\nThe keyword arguments SNR is not expressed in decibels, but as the expected SNR variance ratio. It must be a positive number.\n\nA slightly different version of this model for generating positive definite  matrices has been proposed in (Congedo et al., 2017b)[🎓];  in the model of Eq. [1]\n\nUΛ_iU^H is the signal term, where the signal is supposed sharing the same coordinates for all matrices,\nφ(V_iΔ_iV_i^H) is a structured noise term, which is different for all matrices\nηI is a white noise term, with same variance for all matrices.\n\nSee also: the aforementioned paper and randΛ (randEigvalsMat).\n\nExamples\n\nusing PosDefManifold\nR=randP(10, df=10, eigvalsSNR=1000) # 1 SDP Matrix of size 10x10 #(1)\nH=randP(ComplexF64, 5, eigvalsSNR=10) # 1 Hermitian Matrix of size 5x5 # (2)\nℛ=randP(10, 1000, eigvalsSNR=100) # 1000 SPD Matrices of size 10x10 # (3)\nusing Plots\nheatmap(Matrix(ℛ[1]), yflip=true, c=:bluesreds)\nℋ=randP(ComplexF64, 20, 1000) # 1000 Hermitian Matrices of size 20x20 # (4)\n\n\n\n\n\n"
},

{
    "location": "signalProcessing/#PosDefManifold.regularize!",
    "page": "signalProcessing.jl",
    "title": "PosDefManifold.regularize!",
    "category": "function",
    "text": "(1) regularize!(P::ℍ; <SNR=10e3>)\n(2) regularize!(𝐏::ℍVector; <SNR=10e3>)\n\nAdd white noise to either\n\n(1) a positive definite matrix P of size nn, or\n(2) a 1d array 𝐏 of k positive definite matrices of size nn, of ℍVector type.\n\nThe added noise improves the matrix conditioning with respect to inversion.  This is used to avoid numerical errors when decomposing these matrices  or when evaluating some functions of their eigevalues such as the log.\n\nA constant value is added to all diagonal elements of (1) P  or (2) af all matrices in 𝐏,  that is, on output:\n\ntextrm(1)hspace2ptPleftarrow P+ηI\n\ntextrm(2)hspace2pt𝐏_ileftarrow 𝐏_i+ηI hspace2pttextrmforhspace2pt i=1k\n\nThe amount of added noise η is determined by the SNR  <keyword argument>, which by default is 10000. This is  such that\n\ntextrm(1)hspace2ptSNR=fracdisplaystyletextrmtr(P)displaystyletextrmtr(ηI)\n\ntextrm(2)hspace2ptSNR=fracdisplaystylesum_i=1^ktextrmtr(𝐏_i)displaystyle khspace1pttextrmtr(ηI)\n\nP in (1) must be flagged as Hermitian. See typecasting matrices.\n\nnote: Nota Bene\nThe keyword argument SNR expresses a SNR (signal-to-noise ratio), and is not expressed in decibels,  but as the SNR variance ratio. It must be a positive number. Differently from function randΛrandEigvalsMat, randλrandEigvals and randPrandPosDefMat, the SNR here is not the expected SNR, but the actual SNR.\n\nSee also: randP (randPosDefMat).\n\nExamples\n\n# (1)\nusing LinearAlgebra, Plots, PosDefManifold\nn=3\nU=randU(n)\n# in Q we will write two matrices, the unregularized and regularized matrix side by side\nQ=Matrix{Float64}(undef, n, n*2)\nP=ℍ(U*Diagonal(randn(n).^2)*U\') # generate a real 3x3 positive matrix\nfor i=1:n, j=1:n Q[i, j]=P[i, j] end\nregularize!(P, SNR=5)\nfor i=1:n, j=1:n Q[i, j+n]=P[i, j] end # the regularized matrix is on the right\nheatmap(Matrix(Q), yflip=true, c=:bluesreds)\n\n# (2)\n𝐏=[ℍ(U*Diagonal(randn(3).^2)*U\') for i=1:5] # 5 real 3x3 positive matrices\nregularize!(𝐏, SNR=1000)\n\nRun a test\n\nusing LinearAlgebra\n𝐏=randP(10, 100, SNR=1000); # 100 real Hermitian matrices\nsignalVar=sum(tr(P) for P in 𝐏);\nregularize!(𝐏, SNR=1000);\nsignalPlusNoiseVar=sum(tr(P) for P in 𝐏);\noutput_snr=signalVar/(signalPlusNoiseVar-signalVar)\n# output_snr should be approx. equal to 1000\n\n\n\n\n\n"
},

{
    "location": "signalProcessing/#PosDefManifold.gram",
    "page": "signalProcessing.jl",
    "title": "PosDefManifold.gram",
    "category": "function",
    "text": "gram(X::Matrix{T}) where T<:RealOrComplex\n\nGiven a generic data matrix X, comprised of real or complex elements,  return the normalized Gram matrix, that is,  the covariance matrix of X  corrected by sample size, but without subtracting the mean.\n\nThe result is flagged as Hermitian.  See typecasting matrices.\n\nnote: Nota Bene\nIf X is wide or square (r<=c) return XX^Hc. If X is tall (r>c)            return X^HXr.\n\nExamples\n\nusing PosDefManifold\nX=randn(5, 150);\nG=gram(X) # => G=X*X\'/150\nX=randn(100, 2);\nF=gram(X); # => G=X\'*X/100\n\n\n\n\n\n"
},

{
    "location": "signalProcessing/#PosDefManifold.trade",
    "page": "signalProcessing.jl",
    "title": "PosDefManifold.trade",
    "category": "function",
    "text": "trade(P::ℍ{T}) where T<:RealOrComplex\n\nGiven a positive definite matrix P, return as a 2-tuple the  trace and the determinant of P.  This is used to plot positive matrices in two dimensions  (TraDe plots: log(trace/n) vs. log(determinant), see exemple here below).\n\nP must be flagged by julia as Hermitian.   See typecasting matrices.\n\nExamples\n\nusing PosDefManifold\nP=randP(3)\nt, d=trade(P)  # equivalent to (t, d)=trade(P)\n\n# TraDe plot\nusing Plots\nk=100\nn=10\n𝐏=randP(n, k, SNR=1000); # 100 real Hermitian matrices\nx=Vector{Float64}(undef, k)\ny=Vector{Float64}(undef, k)\nfor i=1:k\n    x[i], y[i] = trade(𝐏[i])\nend\nx=log.(x./n)\ny=log.(y)\nplot(x, y, seriestype=:scatter)\n\n\n\n\n\n"
},

{
    "location": "signalProcessing/#signalProcessing.jl-1",
    "page": "signalProcessing.jl",
    "title": "signalProcessing.jl",
    "category": "section",
    "text": "This unit contains miscellaneous signal processing functions useful in relation to the Riemannian geometry of the manifold of Symmetric Positive Definite (SPD) or Hermitian Positive Definite (HPD) matrices. In Julia those are Hermitian matrices, see typecasting matrices.Function Description\nrandChi², randχ² Generate a random variable distributed as a chi-squared\nrandEigvals, randλ Generate a random vectors of real positive eigenvalues\nrandEigvalsMat, randΛ Generate a random diagonal matrix of real positive eigenvalues\nrandUnitaryMat, randU Generate a random orthogonal or unitary matrix\nrandPosDefMat, randP Generate one or an array of random positive definite matrices\nregularize! Regularize an array of positive definite matrices\ngram Gram matrix of a matrix\ntrade trace and determinant of a matrix as a 2-tuple⋅randChi²\nrandEigvals\nrandEigvalsMat\nrandUnitaryMat\nrandPosDefMat\nregularize!\ngram\ntrade"
},

{
    "location": "classification/#",
    "page": "classification.jl",
    "title": "classification.jl",
    "category": "page",
    "text": ""
},

{
    "location": "classification/#classification.jl-1",
    "page": "classification.jl",
    "title": "classification.jl",
    "category": "section",
    "text": "Unit for classification and related functions.Category Output\n1. Probability functions taking probabilities as input or output"
},

{
    "location": "classification/#PosDefManifold.softmax",
    "page": "classification.jl",
    "title": "PosDefManifold.softmax",
    "category": "function",
    "text": "softmax(χ::Vector{T}) where T<:Real\n\nGiven a real vector of k non-negative scores χ=c_1c_k,  return the vector π=p_1p_k of their  softmax probabilities,  as per\n\np_i=fractextrme^c_isum_i=1^ktextrme^c_i.\n\nExamples\n\nχ=[1.0, 2.3, 0.4, 5.0]\nπ=softmax(χ)\n\n\n\n\n\n"
},

{
    "location": "classification/#Probability-1",
    "page": "classification.jl",
    "title": "Probability",
    "category": "section",
    "text": "Function Description\nsoftmax compute softmax probabilitiesDocumentation:softmax"
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
    "text": "Most functions in PosDefManifold are tested, both for real and complex data input. This unit declares the function testall() that performs all tests.Some functions are fully tested, the others are just executed. Unce you ran it, for each method of each function, a ⭐ sign is printed if the test is succesful, while a ⛔ sign is printed if the test is not succesful. A ☆ sign is printed if the function has been executed correctly.Tests on functions for which a multi-threated version exist are indicated by symbol ( ⏩ ).If there are fails, the concerned functions will be listed as warnings.Note that the first time you execute the test it will take some time as the code will be compiled.This here below is the output of the testall() function (v0.1.3) run on the 20th of May 2019:⭐ PosDefManifold testing utility⭐Starting tests...Unit \'linearAlgebra.jl\'typeofMatrix:           ☆ dim:                    ☆ det1:                   ⭐ ⭐ function tr1:           ⭐ ⭐ normalizeCol!:          ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ispos:                  ⭐ ⭐ colProd:                ⭐ ⭐ ⭐ ⭐ colNorm:                ⭐ ⭐ sumOfSqr:               ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ sumOfSqrDiag:           ⭐ ⭐ ⭐ sumOfSqrTril:           ⭐ ⭐ tr:                     ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ quadraticForm:          ⭐ ⭐ ⭐ ⭐ ⭐ fidelity:               ☆ ☆ fDiag:                  ⭐ DiagOfProd:             ⭐ ⭐ mgs:                    ⭐ ⭐ fVec:                   ⭐ ⭐ ⭐ ⭐ evd:                    ⭐ ⭐ spectralFunctions:      ☆ ☆ ☆ pow:                    ⭐ ⭐ ⭐ invsqrt:                ⭐ ⭐ ⭐ sqr:                    ⭐ ⭐ ⭐ powerIterations:        ⭐ ⭐ ⭐ ⭐ ⭐ choL:                   ⭐ ⭐ ⭐Unit \'signalProcessing.jl\'randλ:                  ☆ randΛ:                  ☆ ☆ randU:                  ⭐ ⭐ randP:                  ☆ ☆ regularize!:            ⭐ ⭐ ⭐ ⭐ gram:                   ☆ ☆ trade:                  ☆ ☆Unit \'riemannianGeometry.jl\'geodesic:               ☆ ☆ ☆ distanceSqr (I):        ☆ ☆ ☆ ☆ distanceSqr (II):       ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ distanceSqr (III):      ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ distance (I):           ☆ ☆ ☆ ☆ distance (II):          ☆ ☆ distanceSqrMat (I):     ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ distanceSqrMat (I ⏩ ):  ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ distanceSqrMat (II):    ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ distanceSqrMat (II ⏩ ): ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ distanceMat (I):        ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ distanceMat (I ⏩ ):     ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ distanceMat (II):       ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ distanceMat (II ⏩ ):    ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ laplacian:              ☆ laplacianEigenMaps:     ☆ spectralEmbedding:      ☆ mean (I):               ☆ ☆ ☆ ☆ ☆ mean (II):              ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ mean (⏩ ):              ☆ ☆ means:                  ☆ ☆ ☆ means (⏩ ):             ☆ ☆ ☆ generalizedMean:        ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ generalizedMean(⏩ ):    ☆ ☆ ☆ ☆ ☆ ☆ geometricMean:          ☆ ☆ ☆ ☆ ☆ ☆ geometricMean(⏩ ):      ☆ ☆ ☆ ☆ ☆ ☆ logdet0Mean:            ⭐ ⭐ ⭐ ⭐ ⭐ logdet0Mean(⏩ ):        ☆ ☆ ☆ ☆ ☆ ☆ wasMean:                ☆ ☆ ☆ ☆ ☆ ☆ wasMean(⏩ ):            ☆ ☆ ☆ ☆ ☆ ☆ powerMean:              ☆ ☆ ☆ ☆ ☆ ☆ powerMean(⏩ ):          ☆ ☆ ☆ ☆ ☆ ☆ logMap:                 ☆ ☆ expMap:                 ☆ ☆ vecP:                   ☆ ☆ matP:                   ☆ ☆ procrustes:             ☆ ☆Unit \'classification.jl\'softmax:                ⭐[ Info: All tests were succesful!"
},

]}

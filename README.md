Guillaume Le Chenadec
Vaibhav Malpani
Rahul Tewari


COMPUTER GRAPHICS: FINAL PROJECT


This is Java code, it runs like the previous assignment.

All our scenes use the same renderer PhotonRenderer.

Some materials have been added: glass, mirror, diffuse and specular, and can be set in the xml file.

Most of the parameters are hardcoded:

 - There is one boolean (useDirect) in Scene.java to specify if we want to use direct lighting (using direct illuminator) or if we want to use only photon mapping.

 - The rest are in src/photon/PhotonMap.java


	public static final int photonsEmitted = 7000000;
	public static final int causticPhotonsNeeded = 1000000;
	private static final int maxDepth = 5;
	private PhotonBoundingVolume photonKDTree = null;
	private PhotonBoundingVolume cPhotonKDTree = null;
	public static final double QUERY_RADIUS_FRACTION = 1.0 / 100.0;  
	public static final double QUERY_RADIUS_FRACTION_CAUSTIC = 1.0 / 100.0;  
	public static final int KNN = 50;
	public static final int KNN_CAUSTIC = 50;
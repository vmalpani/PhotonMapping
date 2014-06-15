package ray.material;

import ray.brdf.BRDF;
import ray.brdf.Lambertian;
import ray.brdf.Specular;
import ray.misc.Color;
import ray.misc.IntersectionRecord;
import ray.misc.LuminaireSamplingRecord;

/**
 * A homogeneous reflecting material, which has the same BRDF at all locations.
 * 
 * @author srm
 */
public class DiffuseAndSpecular implements Material {
	
	BRDF brdf = new Lambertian();
	BRDF specular = new Specular();
	double refractiveIndex = 0.0;

	
	public DiffuseAndSpecular() { }

	public void setBRDF(BRDF brdf) { this.brdf = brdf; }

	public BRDF getBRDF(IntersectionRecord iRec) {
		return brdf;
	}
	
	public void setSpecular(BRDF brdf) { this.specular = brdf; }

	public void setRefractiveIndex(double index) { this.refractiveIndex = index; }

	public BRDF getSpecular(IntersectionRecord iRec) {
		return specular;
	}

	public void emittedRadiance(LuminaireSamplingRecord lRec, Color outRadiance) {
		outRadiance.set(0, 0, 0);
	}

	public boolean isEmitter() {
		return false;
	}

	public boolean isMirror() {
		// TODO Auto-generated method stub
		return false;
	}

	public double getRefractiveIndex() {
		return refractiveIndex;
	}

}

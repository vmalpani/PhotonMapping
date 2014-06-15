package ray.material;

import ray.brdf.BRDF;
import ray.brdf.Lambertian;
import ray.brdf.Specular;
import ray.misc.Color;
import ray.misc.IntersectionRecord;
import ray.misc.LuminaireSamplingRecord;

/**
 * A mirror reflecting material
 * 
 * @author rt2520
 */
public class Mirror implements Material {
	
	BRDF brdf = new Lambertian();
	BRDF specular = new Specular();

	
	public Mirror() { }

	public void setBRDF(BRDF brdf) { this.brdf = brdf; }

	public BRDF getBRDF(IntersectionRecord iRec) {
		return brdf;
	}

	public void setSpecular(BRDF brdf) { this.specular = brdf; }
	//public void setRefractiveIndex(double index) { this.refractiveIndex = index; }

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
		return true;
	}
	
	public double getRefractiveIndex() {
		return 0.0;
	}

}

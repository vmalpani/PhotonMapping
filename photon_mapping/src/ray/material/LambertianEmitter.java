package ray.material;

import ray.brdf.BRDF;
import ray.brdf.Lambertian;
import ray.brdf.Specular;
import ray.misc.Color;
import ray.misc.IntersectionRecord;
import ray.misc.LuminaireSamplingRecord;

/**
 * A Lambertian emitter, which emits the same radiance at all points and in all directions.
 * 
 * @author srm
 */
public class LambertianEmitter implements Material {
	
	Color radiance = new Color();
	BRDF brdf = new Lambertian(new Color(0, 0, 0));
	BRDF specular = new Specular();

	
	public LambertianEmitter() { }
	
	public void setBRDF(BRDF brdf) { this.brdf = brdf; }
	public void setRadiance(Color emittedRadiance) { this.radiance.set(emittedRadiance); }
//	public void setRefractiveIndex(double index) { this.refractiveIndex = index; }

	public BRDF getBRDF(IntersectionRecord iRec) {
		return brdf;
	}

	public void setSpecular(BRDF brdf) { this.specular = brdf; }

	public BRDF getSpecular(IntersectionRecord iRec) {
		return specular;
	}
	
	public void emittedRadiance(LuminaireSamplingRecord lRec, Color outRadiance) {
		outRadiance.set(radiance);
	}

	public boolean isEmitter() {
		return true;
	}
	
	public double getIntensity()
	{
		//return Math.sqrt(radiance.r * radiance.r + radiance.g * radiance.g + radiance.b * radiance.b);
		return (radiance.r * radiance.r + radiance.g * radiance.g + radiance.b * radiance.b)/3;
	}

	public boolean isMirror() {
		// TODO Auto-generated method stub
		return false;
	}
	
	public double getRefractiveIndex() {
		return 0.0;
	}

}

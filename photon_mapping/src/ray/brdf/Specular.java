package ray.brdf;

import ray.math.Frame3;
import ray.math.Point2;
import ray.math.Vector3;
import ray.math.Geometry;
import ray.misc.Color;

/**
 * A Lambertian (constant) BRDF, which has value R/pi where R is the reflectance.
 * 
 * @author srm
 */
public class Specular implements BRDF {
    
    // The material's diffuse reflectance (the fraction of incident irradiance reflected, for any incident distribution).
    Color reflectance = new Color(0.0, 0.0, 0.0);
    
    // For the benefit of the parser
    public Specular() { }
    public void setReflectance(Color reflectance) { this.reflectance.set(reflectance); }
    
    public Specular(Color reflectance) { this.reflectance.set(reflectance); }

    public void evaluate(Frame3 frame, Vector3 incDir, Vector3 reflDir, Color outBRDFValue) {
    	Vector3 outDir = new Vector3(frame.w);
        outDir.scale(-2.0 * (incDir.dot(frame.w)));
        outDir.add(incDir);
        
        if (reflDir.dot(outDir) > 0.9)
        {
        	outBRDFValue.set(reflectance);
        	double cosThetaI = incDir.dot(frame.w);
        	double sinThetaI = Math.sqrt(1 - cosThetaI * cosThetaI);
        	outBRDFValue.scale(1.0 / (sinThetaI * cosThetaI));
        }
    }
    
    public Color getReflectance()
    {
    	return reflectance;
    }

    public void generate(Frame3 frame, Vector3 fixedDir, Vector3 dir, Point2 seed, Color outWeight) {
        Geometry.squareToPSAHemisphere(seed, dir);
        frame.frameToCanonical(dir);
        outWeight.set(reflectance);
    }

    /**
     * @param frame frame comes from IntersectionRecord instance, where w component of this frame align with
     *        the surface normal. 
     * @see ray.brdf.BRDF#pdf(ray.math.Frame3, ray.math.Vector3, ray.math.Vector3)
     */
    public double pdf(Frame3 frame, Vector3 fixedDir, Vector3 dir) {
        return fixedDir.dot(frame.w) / Math.PI;
    }
}

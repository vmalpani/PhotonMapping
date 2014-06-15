package ray.renderer;

import java.util.ArrayList;

import ray.brdf.BRDF;
import ray.material.Material;
import ray.math.Frame3;
import ray.math.Geometry;
import ray.math.Point2;
import ray.math.Point3;
import ray.math.Vector3;
import ray.misc.Color;
import ray.misc.IntersectionRecord;
import ray.misc.LuminaireSamplingRecord;
import ray.misc.Ray;
import ray.misc.Scene;
import ray.surface.Sphere;
import ray.surface.Surface;


/**
 * This class computes direct illumination at a surface by the simplest possible approach: it estimates
 * the integral of incident direct radiance using Monte Carlo integration with a uniform sampling
 * distribution.
 * 
 * The class has two purposes: it is an example to serve as a starting point for other methods, and it
 * is a useful base class because it contains the generally useful <incidentRadiance> function.
 * 
 * @author srm, Changxi Zheng(+)
 */
public class ProjSolidAngleIlluminator extends DirectIlluminator {
    
    
    public void directIllumination(Scene scene, Vector3 incDir, Vector3 outDir, 
            IntersectionRecord iRec, Point2 seed, Color outColor) {
        // W4160 TODO (A)
    	// This method computes a Monte Carlo estimate of reflected radiance due to direct illumination.  It
        // generates samples uniformly wrt. the projected solid angle measure:
        //
        //    f = brdf * radiance
        //    p = 1 / pi
        //    g = f / p = brdf * radiance * pi
        //
        // The same code could be interpreted as an integration wrt. solid angle, as follows:
        //
        //    f = brdf * radiance * cos_theta
        //    p = cos_theta / pi
        //    g = f / p = brdf * radiance * pi
    	// 
    	// As a hint, here are a few steps when I code this function
    	// 1. Generate a random incident direction according to proj solid angle
        //    pdf is constant 1/pi
    	// 2. Find incident radiance from that direction
    	// 3. Estimate reflected radiance using brdf * radiance / pdf = pi * brdf * radiance
    	
    	
    	//Generate a random incident direction incDir and convert it to frame coordinates
        Geometry.squareToPSAHemisphere(seed, incDir);
        iRec.frame.frameToCanonical(incDir);
        //iRec.frame.canonicalToFrame(incDir);
        iRec.frame.canonicalToFrame(outDir);
        
        //Get the brdf of the material where the ray(from the camera) intersects it
        Material material = iRec.surface.getMaterial();
        BRDF brdf = material.getBRDF(iRec);
        Color brdfVal = new Color();
        if ( brdf != null ) {
            brdf.evaluate(iRec.frame, incDir, outDir, brdfVal);
	        // compute incident radiance and scale by the inverse of pdf and reflectance(brdf)
	        scene.incidentRadiance(iRec.frame.o, incDir, outColor);
	        outColor.scale(brdfVal);
	        outColor.scale(Math.PI);
	        return;
        }
        //outColor.set(0, 0, 0);
    }
    
    /**
     * Trace a ray to find direct radiance incident from a particular direction.
     *
     * @param o    the point on surface, where the ray light will arrive
     * @param dir  The direction from which to look for radiance (surface coordinates)
     * @param outRadiance  The radiance found
     */
    public void incidentRadiance(Scene scene, Point3 o, Vector3 dir, Color outRadiance, IntersectionRecord lIntRec) {
        // Trace a ray to find incident (direct) radiance
        Ray ray = new Ray(o, dir);
        ray.makeOffsetRay();
        
        Material material = null;
        if ( scene.getFirstIntersection(lIntRec, ray) && 
            (material = lIntRec.surface.getMaterial()).isEmitter() ) {
            // Hit something -- ask it what its emitted radiance is in our direction
            LuminaireSamplingRecord lSampRec = new LuminaireSamplingRecord();
            lSampRec.set(lIntRec);
            lSampRec.emitDir.set(ray.direction);
            lSampRec.emitDir.scale(-1);
            material.emittedRadiance(lSampRec, outRadiance);
            
            return;
        } 
        // Hit nothing -- background is not direct illumination so return zero
        outRadiance.set(0, 0, 0);
    }
}

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


public class SphereIlluminator extends DirectIlluminator {
    
    
    public void directIllumination(Scene scene, Vector3 incDir, Vector3 outDir, 
            IntersectionRecord iRec, Point2 seed, Color outColor) {
        
    	//Get Spherical lights
    	ArrayList<Surface> lights = scene.getSphericalLights();
    	for (Surface surface : lights)
    	{
    		Sphere lightSphere = (Sphere) surface;
    		Color lightColor = new Color();
	    	//Get the vector from the intersection point to the center of the spherical light source
	    	Point3 center = new Point3();
	    	lightSphere.getCenter(center);
	    	Vector3 pointToCenter = new Vector3();
	    	pointToCenter.sub(center, iRec.frame.o);
	    	
	    	//Generate a frame with w pointing along the point to center vector
	    	//and v pointing along w cross normal
	    	Frame3 samplingFrame = new Frame3();
	    	Vector3 uAxis = new Vector3();
	    	uAxis.cross(pointToCenter, iRec.frame.w);
	    	samplingFrame.set(iRec.frame.o, uAxis, new Vector3(0, 1, 0), pointToCenter);
	    	samplingFrame.initFromWU();
	    	
	    	double lightRadius = lightSphere.getRadius();
	    	double cosAlphaMax = Math.sqrt(1 - (lightRadius * lightRadius) / pointToCenter.squaredLength());
	    	double cosAlpha = 1 - seed.x + seed.x * cosAlphaMax;
	    	double sinAlpha = Math.sqrt(1 - cosAlpha * cosAlpha);
	    	double phi = 2 * Math.PI * seed.y;
	        
	        //Get the brdf of the material where the ray(from the camera) intersects it
	        Material material = iRec.surface.getMaterial();
	        incDir.set(Math.cos(phi) * sinAlpha, Math.sin(phi) * sinAlpha, cosAlpha);
	    	samplingFrame.frameToCanonical(incDir);
	    	//iRec.frame.frameToCanonical(incDir); 
	    	if (incDir.dot(iRec.frame.w) > 0)
	    	{
		    	
		    	iRec.frame.canonicalToFrame(outDir);
	        	BRDF brdf = material.getBRDF(iRec);
	        	Color brdfVal = new Color();
		        if ( brdf != null ) {
		            brdf.evaluate(iRec.frame, incDir, outDir, brdfVal);
			        // compute incident radiance and scale by the inverse of pdf and reflectance(brdf)
		            IntersectionRecord lIntRec = new IntersectionRecord();
			        this.incidentRadiance(scene, lightSphere, iRec.frame.o, incDir, lightColor, lIntRec);
			        lightColor.scale(brdfVal);
			        //Scale with the reciprocal of the pdf of the sampling distribution
			        Vector3 intPointToluminPoint = new Vector3();
			        intPointToluminPoint.sub(lIntRec.frame.o, iRec.frame.o);
			        double oneOverPdf = 2 * Math.PI * intPointToluminPoint.squaredLength() * (1 - cosAlphaMax);
			        lightColor.scale(oneOverPdf);
			        outColor.add(lightColor);		        
		        }
	    	}
    	}
    }
    
    /**
     * Trace a ray to find direct radiance incident from a particular direction.
     *
     * @param o    the point on surface, where the ray light will arrive
     * @param dir  The direction from which to look for radiance (surface coordinates)
     * @param outRadiance  The radiance found
     */
    public void incidentRadiance(Scene scene, Sphere light, Point3 o, Vector3 dir, Color outRadiance, IntersectionRecord lIntRec) {
        // Trace a ray to find incident (direct) radiance
        Ray ray = new Ray(o, dir);
        ray.makeOffsetRay();
        boolean intersectsLight = light.intersect(lIntRec, ray);
        
        IntersectionRecord sceneIntRec = new IntersectionRecord();
        scene.getFirstIntersection(sceneIntRec, ray);
        
        if (intersectsLight  && Math.abs(sceneIntRec.t - lIntRec.t) < Ray.EPSILON) {
            // Hit the luminaire we were trying to sample
        	Material material = lIntRec.surface.getMaterial();
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

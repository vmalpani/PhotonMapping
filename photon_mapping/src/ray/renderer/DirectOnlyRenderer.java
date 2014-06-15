package ray.renderer;

import ray.material.Material;
import ray.math.Point2;
import ray.math.Vector3;
import ray.misc.Color;
import ray.misc.IntersectionRecord;
import ray.misc.LuminaireSamplingRecord;
import ray.misc.Ray;
import ray.misc.Scene;
import ray.sampling.SampleGenerator;

/**
 * A renderer that computes radiance due to emitted and directly reflected light only.
 * 
 * @author srm
 */

public class DirectOnlyRenderer implements Renderer {
    
    /**
     * This is the object that is responsible for computing direct illumination.
     */
    DirectIlluminator direct = null;
        
    /**
     * The default is to compute using uninformed sampling wrt. projected solid angle over the hemisphere.
     */
    public DirectOnlyRenderer() {
        this.direct = new ProjSolidAngleIlluminator();
    }
    
    
    /**
     * This allows the rendering algorithm to be selected from the input file by substituting an instance
     * of a different class of DirectIlluminator.
     * @param direct  the object that will be used to compute direct illumination
     */
    public void setDirectIlluminator(DirectIlluminator direct) {
        this.direct = direct;
    }

    
    public void rayRadiance(Scene scene, Ray ray, SampleGenerator sampler, int sampleIndex, Color outColor) {
        // W4160 TODO (A)
    	// In this function, you need to implement your direct illumination rendering algorithm
    	//
    	// you need:
    	// 1) compute the emitted light radiance from the current surface if the surface is a light surface
    	// 2) direct reflected radiance from other lights. This is by implementing the function
    	//    ProjSolidAngleIlluminator.directIlluminaiton(...), and call direct.directIllumination(...) in this
    	//    function here.
    	// find if the ray intersect with any surface
		IntersectionRecord iRec = new IntersectionRecord();
		Vector3 outDir = new Vector3(ray.direction);
		outDir.scale(-1);
		outDir.normalize();
		outColor.set(0, 0, 0);
		Point2 directSeed = new Point2();
		if (scene.getFirstIntersection(iRec, ray)) {
			if (iRec.surface.getMaterial().isEmitter())
			{
				this.emittedRadiance(iRec, outDir, outColor);
			}
			else //if (!iRec.surface.getMaterial().isMirror())
			{
				for(int sampleIndex2 = 0; sampleIndex2 < sampler.getNumSamples(); sampleIndex2++)
				{
					Color rayColor = new Color();
                    sampler.sample(1, sampleIndex2, directSeed);
                    direct.directIllumination(scene, outDir, iRec, directSeed, rayColor);
                    outColor.add(rayColor);
				}
				outColor.scale(1.0 / sampler.getNumSamples());
			}
			if (iRec.surface.getMaterial().isMirror())
			{
				//Send a ray in the reflected direction
				Vector3 view = new Vector3(outDir);
	        	double viewDotn = iRec.frame.w.dot(view);
	        	Vector3 incDir = new Vector3(iRec.frame.w);
	        	incDir.scale(2 * viewDotn);
	        	incDir.sub(view);
	        	incDir.normalize();

	        	Ray reflectedRay = new Ray(iRec.frame.o, incDir);
	        	reflectedRay.makeOffsetRay();
	        	IntersectionRecord reflectedIntersection = new IntersectionRecord();
	        	if (scene.getFirstIntersection(reflectedIntersection, reflectedRay))
	        	{
		        	Color mirrorColor = new Color();
		        	rayRadiance(scene, reflectedRay, sampler, sampleIndex, mirrorColor);
		        	outColor.scale(mirrorColor);
	        	}
			}
            return;
		}
		scene.getBackground().evaluate(ray.direction, outColor);
    }

    
    /**
     * Compute the radiance emitted by a surface.
     * @param iRec      Information about the surface point being shaded
     * @param dir          The exitant direction (surface coordinates)
     * @param outColor  The emitted radiance is written to this color
     */
    protected void emittedRadiance(IntersectionRecord iRec, Vector3 dir, Color outColor) {
    	// W4160 TODO (A)
        // If material is emitting, query it for emission in the relevant direction.
        // If not, the emission is zero.
    	// This function should be called in the rayRadiance(...) method above
    	LuminaireSamplingRecord lRec = new LuminaireSamplingRecord();
		lRec.emitDir.set(dir);
		lRec.frame.set(iRec.frame);
		lRec.surface = iRec.surface;
		iRec.surface.getMaterial().emittedRadiance(lRec, outColor);
    }
}

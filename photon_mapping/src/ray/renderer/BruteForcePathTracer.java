package ray.renderer;

import ray.brdf.BRDF;
import ray.material.Material;
import ray.math.Geometry;
import ray.math.Point2;
import ray.math.Vector3;
import ray.misc.Color;
import ray.misc.IntersectionRecord;
import ray.misc.LuminaireSamplingRecord;
import ray.misc.Ray;
import ray.misc.Scene;
import ray.sampling.SampleGenerator;

public class BruteForcePathTracer extends PathTracer {
    /**
     * @param scene
     * @param ray
     * @param sampler
     * @param sampleIndex
     * @param outColor
     */
    protected void rayRadianceRecursive(Scene scene, Ray ray, 
            SampleGenerator sampler, int sampleIndex, int level, Color outColor) {
    	// W4160 TODO (B)
    	//
        // Find the visible surface along the ray, then add emitted and reflected radiance
        // to get the resulting color.
    	//
    	// If the ray depth is less than the limit (depthLimit), you need
    	// 1) compute the emitted light radiance from the current surface if the surface is a light surface
    	// 2) reflected radiance from other lights and objects. You need recursively compute the radiance
    	//    hint: You need to call gatherIllumination(...) method.
    	IntersectionRecord iRec = new IntersectionRecord();
    	Vector3 outDir = new Vector3(ray.direction);
		outDir.scale(-1);
		outColor.set(0.0, 0.0, 0.0);
    	if (scene.getFirstIntersection(iRec, ray))
    	{
    		if (iRec.surface.getMaterial().isEmitter())
    		{
    			Color directColor = new Color();
    			this.emittedRadiance(iRec, outDir, directColor);
				outColor.add(directColor);
    		}
    		if (level < this.depthLimit)
	    	{
	    		Color indirectColor = new Color();
	    		gatherIllumination(scene, outDir, iRec, sampler, sampleIndex, level, indirectColor);
	    		outColor.add(indirectColor);
	    	}
    	}
    	else
    		scene.getBackground().evaluate(ray.direction, outColor);    	
    }
}

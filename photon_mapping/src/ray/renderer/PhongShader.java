package ray.renderer;

import java.util.ArrayList;

import ray.brdf.BRDF;
import ray.light.PointLight;
import ray.material.Material;
import ray.math.Geometry;
import ray.math.Point2;
import ray.math.Point3;
import ray.math.Vector3;
import ray.misc.Color;
import ray.misc.IntersectionRecord;
import ray.misc.Ray;
import ray.misc.Scene;
import ray.sampling.SampleGenerator;

public class PhongShader implements Renderer {
    
    private double phongCoeff = 1.5;
    
    public PhongShader() { }
    
    public void setAlpha(double a) {
        phongCoeff = a;
    }
    
    public void rayRadiance(Scene scene, Ray ray, SampleGenerator sampler,
            int sampleIndex, Color outColor) {
        // W4160 TODO (A)
        // Here you need to implement the basic phong reflection model to calculate
        // the color value (radiance) along the given ray. The output color value 
        // is stored in outColor. 
        // 
        // For such a simple rendering algorithm, you might not need Monte Carlo integration
        // In this case, you can ignore the input variable, sampler and sampleIndex.
    	
    	// We have to calculate the dot product of half vector and normal vector
    	// The eye vector here is the input 'ray' reversed, the vector to the light can be 
    	// calculated assuming only point light sources. The intersection record returned by the 
    	// intersect() method gives us the normal at the point of intersection
    	
		// find if the ray intersects with any surface
    	outColor.set(0, 0, 0);
		IntersectionRecord iRec = new IntersectionRecord();
		if (scene.getFirstIntersection(iRec, ray)) {
			ArrayList<PointLight> pointLights = scene.getPointLights();
			for (PointLight pointLight : pointLights)
			{
				// The variable names used here represent the names used in the textbook
				 
				Material intersectionMaterial = iRec.surface.getMaterial();
		    	
		    	Vector3 n = new Vector3(iRec.frame.w);
		    	
		    	Vector3 l = new Vector3();
		    	l.sub(pointLight.location, iRec.frame.o);
		    	l.normalize();
		    	
		    	//Eye vector
		    	Vector3 e = new Vector3(ray.direction);
				e.scale(-1.0);
		    	
		    	double nDotl = n.dot(l);			
		    	if (nDotl > 0)
		    	{
		    		Color diffColor = new Color();
		    		intersectionMaterial.getBRDF(iRec).evaluate(iRec.frame, l, e, diffColor);
		    		diffColor.scale(nDotl);
		    		diffColor.scale(pointLight.diffuse);
		    		Color ambientColor = new Color();
		    		//Homogeneous background, direction is irrelevant 
		    		scene.getBackground().evaluate(new Vector3(), ambientColor);
		    		diffColor.add(ambientColor);
			    	outColor.add(diffColor);
			    	
			    	// Half vector
			    	iRec.frame.canonicalToFrame(e);
			    	iRec.frame.canonicalToFrame(l);
			    	iRec.frame.canonicalToFrame(n);
			    	Vector3 h = new Vector3(e);
			    	h.add(l);
			    	h.normalize();
			    	
			    	Color specColor = new Color();
			    	// TODO: Replace constant specular coeff with material.brdf.specular.b
			    	// Replace 50 with material.shininess
			    	double nDothPowerp  = Math.pow(n.dot(h) , this.phongCoeff);
			    	specColor.b = pointLight.specular.b * 0.5 * nDothPowerp;
			    	specColor.g = pointLight.specular.g * 0.5 * nDothPowerp;
			    	specColor.r = pointLight.specular.r * 0.5 * nDothPowerp;
			    	outColor.add(specColor);
		    	}
			}
			outColor.scale(1.0 / pointLights.size());	    	            
            return;
		}
		scene.getBackground().evaluate(ray.direction, outColor);
    	
    }
}

package photon;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.PriorityQueue;

import photon.accel.PhotonBoundingVolume;
import photon.accel.PhotonBoundingVolume.HeapNode;

import ray.accel.AxisAlignedBoundingBox;
import ray.light.PointLight;
import ray.material.LambertianEmitter;
import ray.math.Geometry;
import ray.math.Point2;
import ray.math.Point3;
import ray.math.Vector3;
import ray.misc.Color;
import ray.misc.IntersectionRecord;
import ray.misc.LuminaireSamplingRecord;
import ray.misc.Ray;
import ray.misc.Scene;
import ray.surface.Surface;
import ray.sampling.IndependentSampler;
import ray.sampling.SampleGenerator;

public class PhotonMap {
	
	
	public static final int photonsEmitted = 5000000;
	public static final int causticPhotonsNeeded = 1000000;
	private static final int maxDepth = 5;
	private PhotonBoundingVolume photonKDTree = null;
	private PhotonBoundingVolume cPhotonKDTree = null;
	public static final double QUERY_RADIUS_FRACTION = 1.0 / 100.0;  
	public static final double QUERY_RADIUS_FRACTION_CAUSTIC = 1.0 / 100.0;  
	public static final int KNN = 50;
	public static final int KNN_CAUSTIC = 50;
	
	
	public void buildPhotonMap(Scene scene)
	{
		ArrayList<Photon> photonMap = new ArrayList<Photon>();
		IndependentSampler sampler = new IndependentSampler();
		Point2 seed = new Point2();
		
		ArrayList<Surface> luminaires = scene.getLuminairesAsSurfaces();
		ArrayList<PointLight> pointLights = scene.getPointLights();
		double totalIntensity = 0;
		for (Surface luminaire : luminaires)
		{
			totalIntensity += ((LambertianEmitter)(luminaire.getMaterial())).getIntensity();
		}
		
		//send photons from the area lights
		for (PointLight pointLight : pointLights)
		{
			totalIntensity += pointLight.intensity;
			
		}
		//send photons from the point light sources
		for (Surface luminaire : luminaires)
		{
	    	double intensity = ((LambertianEmitter)(luminaire.getMaterial())).getIntensity();
			for (int i = 0; i < photonsEmitted * (intensity / totalIntensity); i++)
			{
				sampler.sample(0,0,seed);
				Point3 dummy = new Point3();
				LuminaireSamplingRecord lRec = new LuminaireSamplingRecord();
				luminaire.chooseSamplePoint(dummy, seed, lRec);
				
				sampler.sample(0, 0, seed);
				Vector3 randomDir = new Vector3();
				Geometry.squareToCosineWeightedHemisphere(seed, randomDir);
				lRec.frame.frameToCanonical(randomDir);
				
				//System.out.println(lRec.frame.w.x + ", " + lRec.frame.w.y + ", " + lRec.frame.w.z + ", ");
				Ray ray = new Ray(lRec.frame.o, randomDir);
				ray.makeOffsetRay();

				Photon photon = new Photon();
				((LambertianEmitter)(luminaire.getMaterial())).emittedRadiance(lRec, photon.color);
				//photon.color.scale(intensity);
				//photon.color.scale(totalIntensity / (photonsEmitted * intensity));
				photon.color.scale(intensity / photonsEmitted);
				tracePhoton(scene, sampler, ray, photon, 0, photonMap);
			}
		}
		for (PointLight pointLight : pointLights)
		{
			for (int i = 0; i < photonsEmitted * (pointLight.intensity / totalIntensity); i++)
			{
				sampler.sample(0,0,seed);
				Vector3 randomDir = new Vector3();
				Geometry.squareToSphere(seed, randomDir);
				
				
				Ray ray = new Ray(pointLight.location, randomDir);
				ray.makeOffsetRay();
				
				/*System.out.println(ray.direction.x);
				System.out.println(ray.direction.y);
				System.out.println(ray.direction.z);
*/
				
				Photon photon = new Photon();
				photon.color.set(pointLight.diffuse);
				//photon.color.scale(pointLight.intensity);
				photon.color.scale(totalIntensity / photonsEmitted);
				tracePhoton(scene, sampler, ray, photon, 0 ,photonMap);
			}
		}
		photonKDTree = new PhotonBoundingVolume(photonMap);
		
	}
	
	public void buildCausticPhotonMap(Scene scene)
	{
		ArrayList<Photon> cPhotonMap = new ArrayList<Photon>();
		IndependentSampler sampler = new IndependentSampler();
		Point2 seed = new Point2();
		
		ArrayList<Surface> luminaires = scene.getLuminairesAsSurfaces();
		ArrayList<PointLight> pointLights = scene.getPointLights();
		double totalIntensity = 0;
		for (Surface luminaire : luminaires)
		{
			totalIntensity += ((LambertianEmitter)(luminaire.getMaterial())).getIntensity();
		}
		
		//send photons from the area lights
		for (PointLight pointLight : pointLights)
		{
			totalIntensity += pointLight.intensity;
			
		}
		
		ArrayList<Surface> specularSurfaces = scene.getSpecularSurfaces();
		
		for (Surface luminaire : luminaires)
		{
			
			int startIndexInPhotonMap = cPhotonMap.size();
			
			int totalNbDirectionsGenerated = 0;
	    	double intensity = ((LambertianEmitter)(luminaire.getMaterial())).getIntensity();
			for (int i = 0; i < causticPhotonsNeeded * (intensity / totalIntensity); i++)
			{
				sampler.sample(0,0,seed);
				Point3 dummy = new Point3();
				LuminaireSamplingRecord lRec = new LuminaireSamplingRecord();
				luminaire.chooseSamplePoint(dummy, seed, lRec);
				
				boolean hitSpecularObject = false;
				while (!hitSpecularObject)
				{
					boolean hitSpecularBB = false;
					Vector3 randomDir = new Vector3();
					Ray ray = new Ray(lRec.frame.o, randomDir);
					while (!hitSpecularBB)
					{
						sampler.sample(0, 0, seed);
						Geometry.squareToCosineWeightedHemisphere(seed, randomDir);
						lRec.frame.frameToCanonical(randomDir);
						
						ray.direction.set(randomDir);
						totalNbDirectionsGenerated++;
						ray.makeOffsetRay();
						
						for (Surface specSurf : specularSurfaces)
						{
							if (specSurf.getBoundingBox().intersect(ray))
							{
								hitSpecularBB = true;
								break;
							}
						}
					}
					Photon photon = new Photon();
					((LambertianEmitter)(luminaire.getMaterial())).emittedRadiance(lRec, photon.color);
					//photon.color.scale(intensity);
					//photon.color.scale(totalIntensity / (causticPhotonsNeeded * intensity));
					photon.color.scale(totalIntensity / causticPhotonsNeeded );
					hitSpecularObject = traceCausticPhoton(scene, sampler, ray, photon, 0, cPhotonMap, false);
				}
			}
			
			//scale the intensity of the photon sent by the total number of photons we "send" (the directions that we generate but that do not hit a specular object)
			double ratio = (causticPhotonsNeeded * (intensity / totalIntensity)) / totalNbDirectionsGenerated;
			
			System.out.println(cPhotonMap.size());
			for(int i = startIndexInPhotonMap; i < startIndexInPhotonMap + (int)(causticPhotonsNeeded * (intensity / totalIntensity)); i++)
				cPhotonMap.get(i).color.scale(ratio);
			
		}
		for (PointLight pointLight : pointLights)
		{
			
			int startIndexInPhotonMap = cPhotonMap.size();
			int totalNbDirectionsGenerated = 0;

			for (int i = 0; i < causticPhotonsNeeded * (pointLight.intensity / totalIntensity); i++)
			{
				boolean hitSpecularObject = false;
				while (!hitSpecularObject)
				{
					Vector3 randomDir = new Vector3();
					boolean hitSpecularBB = false;
					Ray ray = new Ray(pointLight.location, randomDir);
					while (!hitSpecularBB)
					{
						sampler.sample(0,0,seed);
						Geometry.squareToSphere(seed, randomDir);
									
						ray.direction.set(randomDir);
						ray.makeOffsetRay();
						totalNbDirectionsGenerated++;

						
						for (Surface specSurf : specularSurfaces)
						{
							if (specSurf.getBoundingBox().intersect(ray))
							{
								hitSpecularBB = true;
								break;
							}
						}
					}
					
					Photon photon = new Photon();
					photon.color.set(pointLight.diffuse);
					//photon.color.scale(pointLight.intensity);
					photon.color.scale(totalIntensity / causticPhotonsNeeded);
					hitSpecularObject = traceCausticPhoton(scene, sampler, ray, photon, 0 ,cPhotonMap, false);
				}
			}
			//scale the intensity of the photon sent by the total number of photons we "send" (the directions that we generate but that do not hit a specular object)
			double ratio =  (causticPhotonsNeeded * (pointLight.intensity / totalIntensity)) / totalNbDirectionsGenerated;
			
			for(int i = startIndexInPhotonMap; i < startIndexInPhotonMap + (int)(causticPhotonsNeeded * (pointLight.intensity / totalIntensity)); i++)
				cPhotonMap.get(i).color.scale(ratio);
		}
		cPhotonKDTree = new PhotonBoundingVolume(cPhotonMap);
		
	}
	
	private void tracePhoton(Scene scene, SampleGenerator sampler, Ray ray, Photon photon, int depth, ArrayList<Photon> photonMap)
	{
		IntersectionRecord iRec = new IntersectionRecord();
		if(scene.getFirstIntersection(iRec, ray))
		{
			Color brdfDiffuse = new Color();	
			brdfDiffuse.set(iRec.surface.getMaterial().getBRDF(iRec).getReflectance());
			
			//if the material has some diffuse component, we store it
			if ((brdfDiffuse.r > 0.0001 || brdfDiffuse.g > 0.0001 || brdfDiffuse.b > 0.0001) && (depth > 0 || !scene.useDirect()))
			{
				Photon intersectionPhoton = new Photon();
				intersectionPhoton.color.set(photon.color);
				intersectionPhoton.position.set(iRec.frame.o);
				intersectionPhoton.incidentDirection.set(ray.direction);
				if (depth == 0)
					intersectionPhoton.isDirect = true;
				photonMap.add(intersectionPhoton);
			}
			
			if(depth < maxDepth)
			{
				Point2 seed = new Point2();
				//ray.direction.scale(-1);
				
				Vector3 outDir = new Vector3();
/*				Geometry.squareToHemisphere(seed, outDir);
	            iRec.frame.frameToCanonical(outDir);
*/
				
				Color brdfSpecular = new Color();
				brdfSpecular.set(iRec.surface.getMaterial().getSpecular(iRec).getReflectance());

				//the diffuse reflectance and specular reflectance have to sum to less than 1
				Color brdfDiffuseTimesPhotonIntensity = new Color(brdfDiffuse);
				brdfDiffuseTimesPhotonIntensity.scale(photon.color);
				double denominator = (photon.color.r > photon.color.g) ? ((photon.color.r > photon.color.b) ? photon.color.r : photon.color.b) 
						: ((photon.color.g > photon.color.b) ? photon.color.g : photon.color.b);
				double numerator = (brdfDiffuseTimesPhotonIntensity.r > brdfDiffuseTimesPhotonIntensity.g) ? ((brdfDiffuseTimesPhotonIntensity.r > brdfDiffuseTimesPhotonIntensity.b) ? brdfDiffuseTimesPhotonIntensity.r : brdfDiffuseTimesPhotonIntensity.b) 
						: ((brdfDiffuseTimesPhotonIntensity.g > brdfDiffuseTimesPhotonIntensity.b) ? brdfDiffuseTimesPhotonIntensity.g : brdfDiffuseTimesPhotonIntensity.b);
				
				double diffuseProbability = numerator/denominator;
				
				Color brdfSpecularTimesPhotonIntensity = new Color(brdfSpecular);
				brdfSpecularTimesPhotonIntensity.scale(photon.color);
				numerator = (brdfSpecularTimesPhotonIntensity.r > brdfSpecularTimesPhotonIntensity.g) ? ((brdfSpecularTimesPhotonIntensity.r > brdfSpecularTimesPhotonIntensity.b) ? brdfSpecularTimesPhotonIntensity.r : brdfSpecularTimesPhotonIntensity.b) 
						: ((brdfSpecularTimesPhotonIntensity.g > brdfSpecularTimesPhotonIntensity.b) ? brdfSpecularTimesPhotonIntensity.g : brdfSpecularTimesPhotonIntensity.b);
				
				double specularProbability = numerator/denominator;
				
				sampler.sample(0,0,seed);
				
				//diffuse reflection
				if(seed.x < diffuseProbability)
				{
					sampler.sample(0,0,seed);
					Geometry.squareToHemisphere(seed, outDir);
		            iRec.frame.frameToCanonical(outDir);	
		            
		            photon.color.set(brdfDiffuseTimesPhotonIntensity);
		            photon.color.scale(1.0/diffuseProbability);
		            
		            ray.origin.set(iRec.frame.o);
		            ray.direction.set(outDir);
					ray.makeOffsetRay();

		            tracePhoton(scene, sampler, ray, photon, depth+1, photonMap);
				}
				//specular reflection
				else if(seed.x < diffuseProbability + specularProbability)
				{
					outDir.set(iRec.frame.w);
			        outDir.scale(-2.0*(ray.direction).dot(iRec.frame.w));
			        outDir.add(ray.direction);

			        photon.color.set(brdfSpecularTimesPhotonIntensity);
		            photon.color.scale(1.0/specularProbability);
		            
		            ray.origin.set(iRec.frame.o);
		            ray.direction.set(outDir);
					ray.makeOffsetRay();

		            tracePhoton(scene, sampler, ray, photon, depth+1, photonMap);
				}
				
			}
		}	
	}
	
	private boolean traceCausticPhoton(Scene scene, SampleGenerator sampler, Ray ray, Photon photon, int depth, ArrayList<Photon> photonMap, boolean insideDielectric)
	{
		
		
		IntersectionRecord iRec = new IntersectionRecord();
		if(scene.getFirstIntersection(iRec, ray))
		{
			Color brdfSpecular = new Color();
			brdfSpecular.set(iRec.surface.getMaterial().getSpecular(iRec).getReflectance());
			
			//System.out.println("ok1");

			if(brdfSpecular.r < 0.0001 && brdfSpecular.g < 0.0001 && brdfSpecular.b < 0.0001 && depth == 0 && iRec.surface.getMaterial().getRefractiveIndex() < Ray.EPSILON )
			{
				return false;
			}

			//System.out.println("ok");
			Color brdfDiffuse = new Color();	
			brdfDiffuse.set(iRec.surface.getMaterial().getBRDF(iRec).getReflectance());
			
			//if the material has some diffuse component, we store the caustic photon and stop the tracing
			if((brdfDiffuse.r > 0.0001 || brdfDiffuse.g > 0.0001 || brdfDiffuse.b > 0.0001) && (depth > 0))
			{
				Photon intersectionPhoton = new Photon();
				intersectionPhoton.color.set(photon.color);
				intersectionPhoton.position.set(iRec.frame.o);
				intersectionPhoton.incidentDirection.set(ray.direction);
				if (depth == 0)
					intersectionPhoton.isDirect = true;
				photonMap.add(intersectionPhoton);
				
				if(Double.isNaN(photon.color.r))
					System.out.println("r nan");
				if(Double.isNaN(photon.color.g))
					System.out.println("g nan");
				if(Double.isNaN(photon.color.b))
					System.out.println("b nan");
				
				return true;
			}
			
			if(depth < maxDepth)
			{
				Point2 seed = new Point2();
				//ray.direction.scale(-1);
				
				Vector3 outDir = new Vector3();
//				Geometry.squareToHemisphere(seed, outDir);
//	            iRec.frame.frameToCanonical(outDir);

				
				

				//the diffuse reflectance and specular reflectance have to sum to less than 1
				Color brdfDiffuseTimesPhotonIntensity = new Color(brdfDiffuse);
				brdfDiffuseTimesPhotonIntensity.scale(photon.color);
				double denominator = (photon.color.r > photon.color.g) ? ((photon.color.r > photon.color.b) ? photon.color.r : photon.color.b) 
						: ((photon.color.g > photon.color.b) ? photon.color.g : photon.color.b);
				double numerator = (brdfDiffuseTimesPhotonIntensity.r > brdfDiffuseTimesPhotonIntensity.g) ? ((brdfDiffuseTimesPhotonIntensity.r > brdfDiffuseTimesPhotonIntensity.b) ? brdfDiffuseTimesPhotonIntensity.r : brdfDiffuseTimesPhotonIntensity.b) 
						: ((brdfDiffuseTimesPhotonIntensity.g > brdfDiffuseTimesPhotonIntensity.b) ? brdfDiffuseTimesPhotonIntensity.g : brdfDiffuseTimesPhotonIntensity.b);
				
				double diffuseProbability = numerator/denominator;
				
				Color brdfSpecularTimesPhotonIntensity = new Color(brdfSpecular);
				brdfSpecularTimesPhotonIntensity.scale(photon.color);
				numerator = (brdfSpecularTimesPhotonIntensity.r > brdfSpecularTimesPhotonIntensity.g) ? ((brdfSpecularTimesPhotonIntensity.r > brdfSpecularTimesPhotonIntensity.b) ? brdfSpecularTimesPhotonIntensity.r : brdfSpecularTimesPhotonIntensity.b) 
						: ((brdfSpecularTimesPhotonIntensity.g > brdfSpecularTimesPhotonIntensity.b) ? brdfSpecularTimesPhotonIntensity.g : brdfSpecularTimesPhotonIntensity.b);
				
				double specularProbability = numerator/denominator;
				
				sampler.sample(0,0,seed);
				
				//diffuse reflection
				if(seed.x < diffuseProbability && depth == 0)
					return false;
				if(seed.x < diffuseProbability)
				{
					sampler.sample(0,0,seed);
					Geometry.squareToHemisphere(seed, outDir);
		            iRec.frame.frameToCanonical(outDir);	
		            
		            photon.color.set(brdfDiffuseTimesPhotonIntensity);
		            photon.color.scale(1.0/diffuseProbability);
		            
		            ray.origin.set(iRec.frame.o);
		            ray.direction.set(outDir);
					ray.makeOffsetRay();

		            return traceCausticPhoton(scene, sampler, ray, photon, depth+1, photonMap, insideDielectric);
				}
				//specular reflection
				else if(seed.x < diffuseProbability + specularProbability)
				{
					outDir.set(iRec.frame.w);
			        outDir.scale(-2.0*(ray.direction).dot(iRec.frame.w));
			        outDir.add(ray.direction);

				//System.out.println("before : " + photon.color.r + "  "  +  photon.color.b + "  "  +  photon.color.g);
			        photon.color.set(brdfSpecularTimesPhotonIntensity);
		            photon.color.scale(1.0/specularProbability);
				//System.out.println("after : " + photon.color.r + "  "  +  photon.color.b + "  "  +  photon.color.g);
		            
		            ray.origin.set(iRec.frame.o);
		            ray.direction.set(outDir);
					ray.makeOffsetRay();

					//System.out.println(ray.direction.x + ", "+ ray.direction.y  + ", "+ ray.direction.z);
		            return traceCausticPhoton(scene, sampler, ray, photon, depth+1, photonMap, insideDielectric);
				}
				double nT = iRec.surface.getMaterial().getRefractiveIndex();
				if (nT > Ray.EPSILON)//refracted
				{
					//Naming conventions from Shirley
					double n = 1;
					Vector3 normal = new Vector3(iRec.frame.w);
					if (insideDielectric)
					{
						nT = nT + n;
						n = nT - n;
						nT = nT -n;
						normal.scale(-1.0);
					}
					Ray refractedRay = new Ray(ray);
					refractedRay.direction.scale(-1.0);
					double cosTheta = normal.dot(refractedRay.direction);
					double cosSqPhi = 1 - ((n * n) * (1 -(cosTheta * cosTheta))) / (nT * nT);
					if (cosSqPhi > 0) // refraction
					{
						double cosPhi = Math.sqrt(cosSqPhi);
						double schlickReflectivity = (nT - n) / (nT + n);
						schlickReflectivity *= schlickReflectivity;
						double schlickReflectivityTheta = schlickReflectivity +  (1 - schlickReflectivity) * Math.pow(1 - (!insideDielectric ? cosTheta : cosPhi), 5);
						
						//double schlickReflectivityTheta = schlickReflectivity +  (1 - schlickReflectivity) * Math.pow(1 - (!insideDielectric ? cosTheta : cosPhi), 5);
						sampler.sample(0,0,seed);
						seed.x = schlickReflectivity + (1 - schlickReflectivity) * seed.x; 
						if (seed.x < schlickReflectivityTheta)
						{
							//reflection
							outDir.set(normal);
					        outDir.scale(-2.0*(ray.direction).dot(normal));
					        outDir.add(ray.direction);

					        ray.origin.set(iRec.frame.o);
				            ray.direction.set(outDir);
							ray.makeOffsetRay();
							
							if(!insideDielectric)
								depth++;

				            return traceCausticPhoton(scene, sampler, ray, photon, depth, photonMap, insideDielectric);

						}
						else
						{
							Vector3 refractedDir = new Vector3(normal);
							refractedDir.scale(cosTheta);
							refractedDir.add(ray.direction);
							refractedDir.scale(n / nT);
							refractedDir.scaleAdd(-cosPhi, normal);
							refractedRay.set(iRec.frame.o, refractedDir);
							refractedRay.makeOffsetRay();
							
							if(!insideDielectric)
								depth++;
				            return traceCausticPhoton(scene, sampler, refractedRay, photon, depth, photonMap, !insideDielectric);

						}
					}
					else
					{	
						//reflection
						outDir.set(normal);
				        outDir.scale(-2.0*(ray.direction).dot(normal));
				        outDir.add(ray.direction);

				        ray.origin.set(iRec.frame.o);
			            ray.direction.set(outDir);
						ray.makeOffsetRay();
						
						if(!insideDielectric)
							depth++;
			            return traceCausticPhoton(scene, sampler, ray, photon, depth, photonMap, insideDielectric);
					 }
				}
				// else absorbed				
			}
		}

		return false;
		
	}
	
	public void getFluxFromKNearest(Point3 point, double sceneRadius,
															PriorityQueue<PhotonBoundingVolume.HeapNode> globalHeap,
															PriorityQueue<PhotonBoundingVolume.HeapNode> causticHeap)
	{
		PhotonBoundingVolume.PhotonDistanceComparator comp = new PhotonBoundingVolume.PhotonDistanceComparator();
		PriorityQueue<PhotonBoundingVolume.HeapNode> heap = new PriorityQueue<PhotonBoundingVolume.HeapNode>(KNN, comp);
		
		double radius = sceneRadius * QUERY_RADIUS_FRACTION;
		
		this.photonKDTree.getFluxFromKNearestRecursive(point, radius * radius, heap);
		globalHeap.addAll(heap);
		heap = new PriorityQueue<PhotonBoundingVolume.HeapNode>(KNN_CAUSTIC, comp);
		radius = sceneRadius * QUERY_RADIUS_FRACTION_CAUSTIC;
		this.cPhotonKDTree.getFluxFromKNearestRecursive(point, radius * radius, heap);
		causticHeap.addAll(heap);
	}
}

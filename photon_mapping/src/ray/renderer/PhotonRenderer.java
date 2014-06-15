package ray.renderer;

import java.util.Iterator;
import java.util.PriorityQueue;

import photon.PhotonMap;
import photon.accel.PhotonBoundingVolume;
import ray.math.Point2;
import ray.math.Vector3;
import ray.misc.Color;
import ray.misc.IntersectionRecord;
import ray.misc.LuminaireSamplingRecord;
import ray.misc.Ray;
import ray.misc.Scene;
import ray.sampling.SampleGenerator;

public class PhotonRenderer implements Renderer{

	public PhotonMap photonMap = null;
	
	private DirectIlluminator direct = null;
	public PhotonRenderer() {
		this.direct = new LuminairesIlluminator();
	}
	
	/**
     * This allows the rendering algorithm to be selected from the input file by substituting an instance
     * of a different class of DirectIlluminator.
     * @param direct  the object that will be used to compute direct illumination
     */
    public void setDirectIlluminator(DirectIlluminator direct) {
        this.direct = direct;
    }
	
    public void rayRadianceDielectric(Scene scene, Ray ray, SampleGenerator sampler,
			int sampleIndex, Color outColor, boolean insideDielectric, int nbInternalReflections)
    {
    		IntersectionRecord iRec = new IntersectionRecord();
		outColor.set(0.0, 0.0, 0.0);
		if (scene.getFirstIntersection(iRec, ray))
		{	
			Vector3 outDir = new Vector3(ray.direction);
			outDir.scale(-1.0);
			Color emitted = new Color();
			emittedRadiance(iRec, ray.direction, emitted);
			outColor.add(emitted);
			if (!iRec.surface.getMaterial().isEmitter() && scene.useDirect())
			{
				Point2 directSeed = new Point2();
				sampler.sample(1, sampleIndex, directSeed);
				this.direct.directIllumination(scene, outDir, iRec, directSeed, outColor);	
				//outColor.scale(1.0/100000.0);
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
				rayRadianceDielectric(scene, reflectedRay, sampler, 0, mirrorColor, insideDielectric, 0);
		        	//rayRadiance(scene, reflectedRay, sampler, sampleIndex, mirrorColor);
		        	outColor.add(mirrorColor);
	        	}
	        	return;
			}
			
			double nT = iRec.surface.getMaterial().getRefractiveIndex();
			if (nT > Ray.EPSILON && nbInternalReflections < 3)//refracted
			{
				//Naming conventions from Shirley
				double n = 1.0;
				Vector3 normal = new Vector3(iRec.frame.w);
				if (insideDielectric)
				{
					//nT = nT + n;
					//n = nT - n;
					//nT = nT -n;
					n = nT;
					nT=1.0;
					normal.scale(-1.0);
				}
				
				Vector3 reflectedDir = new Vector3(normal);
				reflectedDir.scale(-2.0*(ray.direction).dot(normal));
				reflectedDir.add(ray.direction);
				
				Ray reflectedRay = new Ray(iRec.frame.o, reflectedDir);
		        reflectedRay.makeOffsetRay();
		        
		        Ray refractedRay = new Ray(ray);
				refractedRay.direction.scale(-1.0);
				double cosTheta = normal.dot(refractedRay.direction);
				//System.out.println("cosTheta " + cosTheta);
				double cosSqPhi = 1 - ((n * n) * (1 -(cosTheta * cosTheta))) / (nT * nT);
				//System.out.println("cosSqPhi " + cosSqPhi);

				Color reflectColor = new Color();
				Color refractedColor = new Color();
				refractedColor.set(0.0,0.0,0.0);

				double schlickReflectivityTheta = 1.0;
				
		    		if (cosSqPhi > 0) // refraction
				{

		    		double cosPhi = Math.sqrt(cosSqPhi);
					double schlickReflectivity = (nT - n) / (nT + n);
					schlickReflectivity *= schlickReflectivity;
					schlickReflectivityTheta = schlickReflectivity +  (1 - schlickReflectivity) * Math.pow(1 - (!insideDielectric ? cosTheta : cosPhi), 5);
					//System.out.println(schlickReflectivityTheta);
					Vector3 refractedDir = new Vector3(normal);
					refractedDir.scale(cosTheta);
					refractedDir.add(ray.direction);
					refractedDir.scale(n / nT);
					Vector3 temp = new Vector3(normal);
					temp.scale(-cosPhi);
					refractedDir.add(temp);
					//refractedDir.scaleAdd(-cosPhi, normal);
					refractedRay.set(iRec.frame.o, refractedDir);
					refractedRay.makeOffsetRay();
					
					//Ray outRay = new Ray();
					//traceRefractedRay(scene, nT, refractedRay, 0, outRay);
				
					rayRadianceDielectric(scene, refractedRay, sampler, 0, refractedColor, !insideDielectric, 0);
					//outColor.scaleAdd(Math.max(0.0, 1 - schlickReflectivityTheta), reflectColor);
				}
				
		    	reflectColor.set(0,0,0);
		    	if(insideDielectric)
		    		nbInternalReflections++;
		    	else
		    		nbInternalReflections = 0;
		    	
				rayRadianceDielectric(scene, reflectedRay, sampler, 0, reflectColor, insideDielectric, nbInternalReflections);

				refractedColor.scale(1.0-schlickReflectivityTheta);
				refractedColor.scaleAdd(schlickReflectivityTheta, reflectColor);
				outColor.add(refractedColor);

				return;
		    	//outColor.scaleAdd(Math.max(0.0, schlickReflectivityTheta), reflectColor);
			}
			else if(nT < Ray.EPSILON)
			{
				
			
				//photon color
				Vector3 extents = scene.getBoundingBoxExtents();
				double sceneRadius = (extents.x > extents.y) ? ((extents.x > extents.z) ? extents.x : extents.z) 
										: ((extents.y > extents.z) ? extents.y : extents.z);
				//System.out.println( "Scene radius : " + sceneRadius);
				PhotonBoundingVolume.PhotonDistanceComparator comp = new PhotonBoundingVolume.PhotonDistanceComparator();
				PriorityQueue<PhotonBoundingVolume.HeapNode> globalHeap = new PriorityQueue<PhotonBoundingVolume.HeapNode>(PhotonMap.KNN, comp);
				PriorityQueue<PhotonBoundingVolume.HeapNode> causticHeap = new PriorityQueue<PhotonBoundingVolume.HeapNode>(PhotonMap.KNN, comp);
				this.photonMap.getFluxFromKNearest(iRec.frame.o, sceneRadius, globalHeap, causticHeap);
				Iterator<PhotonBoundingVolume.HeapNode> gIter = globalHeap.iterator();
				Iterator<PhotonBoundingVolume.HeapNode> cIter = causticHeap.iterator();
				//if(heap.size() < 2)
					//System.out.println(heap.size());
				
				//!!reverse incident directions for brdf?
				ray.direction.scale(-1.0);
				
				double maxRadius = 10.0;
				if(globalHeap.size() > 0)
					maxRadius = globalHeap.peek().distanceFromPoint;
				
				while (gIter.hasNext())
				{
					PhotonBoundingVolume.HeapNode heapNode = gIter.next();
					Color photonContribution = new Color(heapNode.photon.color);
					Color brdf = new Color();
					
					Vector3 incDir = new Vector3(heapNode.photon.incidentDirection);
					incDir.scale(-1.0);
					if(incDir.dot(iRec.frame.w) > 0)
					{
						iRec.surface.getMaterial().getBRDF(iRec).evaluate(iRec.frame, incDir, ray.direction, brdf);
	//					if (heapNode.photon.isDirect)
	//					{
							Color specBRDF = new Color();
							iRec.surface.getMaterial().getSpecular(iRec).evaluate(iRec.frame, incDir, ray.direction, specBRDF);
							brdf.add(specBRDF);
	//					}
						photonContribution.scale(brdf);
						//TODO: Each radius or max radius
						
						//the distance stored in distance from point is already squared

						//photonContribution.scale(1 / (heapNode.distanceFromPoint * Math.PI));
						//photonContribution.scale(1 / (heapNode.distanceFromPoint * Math.PI));
						photonContribution.scale(1 / (maxRadius * Math.PI));

						outColor.add(photonContribution);
					}
				}
	
				//ray.direction.scale(-1.0);
				if(causticHeap.size() > 0)
					maxRadius = causticHeap.peek().distanceFromPoint;
				
				while (cIter.hasNext())
				{
					PhotonBoundingVolume.HeapNode heapNode = cIter.next();
					Color photonContribution = new Color(heapNode.photon.color);
					
					
					
					Color brdf = new Color();
					Vector3 incDir = new Vector3(heapNode.photon.incidentDirection);
					incDir.scale(-1.0);
					if(incDir.dot(iRec.frame.w) > 0)
					{
						iRec.surface.getMaterial().getBRDF(iRec).evaluate(iRec.frame, incDir, ray.direction, brdf);
	//					if (heapNode.photon.isDirect)
	//					{
							Color specBRDF = new Color();
							iRec.surface.getMaterial().getSpecular(iRec).evaluate(iRec.frame, incDir, ray.direction, specBRDF);
							brdf.add(specBRDF);
	//					}
						photonContribution.scale(brdf);
						
						//TODO: Each radius or max radius
						
						//the distance stored in distance from point is already squared
						if (heapNode.distanceFromPoint > 0)
							photonContribution.scale(1 / (maxRadius * Math.PI));

					//		photonContribution.scale(1 / (heapNode.distanceFromPoint * Math.PI));

					//		photonContribution.scale(1 / (heapNode.distanceFromPoint * heapNode.distanceFromPoint * Math.PI));
						outColor.add(photonContribution);
						//if(Double.isNaN(photonContribution.r) || Double.isNaN(photonContribution.g) || Double.isNaN(photonContribution.b) )
							//System.out.println("pb");
					}
				}
				return;
				
			}			
			else
				return;	
		}
		
		scene.getBackground().evaluate(ray.direction, outColor);
    
    }
    
    
    
	public void rayRadiance(Scene scene, Ray ray, SampleGenerator sampler, int sampleIndex, Color outColor) 
	{
		
		this.rayRadianceDielectric(scene, ray, sampler, sampleIndex, outColor, false, 0);
		/*
		IntersectionRecord iRec = new IntersectionRecord();
		outColor.set(0.0, 0.0, 0.0);
		if (scene.getFirstIntersection(iRec, ray))
		{		
			Vector3 outDir = new Vector3(ray.direction);
			outDir.scale(-1.0);
			
			if (scene.useDirect())
			{
				Point2 directSeed = new Point2();
				sampler.sample(1, sampleIndex, directSeed);
				this.direct.directIllumination(scene, outDir, iRec, directSeed, outColor);			
			}
			
			double nT = iRec.surface.getMaterial().getRefractiveIndex();
			if (nT > Ray.EPSILON)//refracted
			{
				//Naming conventions from Shirley
				double n = 1;
				Vector3 normal = new Vector3(iRec.frame.w);
				
				Vector3 reflectedDir = new Vector3(normal);
				reflectedDir.scale(-2.0*(ray.direction).dot(normal));
				reflectedDir.add(ray.direction);
				
				Ray reflectedRay = new Ray(iRec.frame.o, reflectedDir);
		        reflectedRay.makeOffsetRay();
		        
		        Ray refractedRay = new Ray(ray);
				refractedRay.direction.scale(-1.0);
				double cosTheta = normal.dot(refractedRay.direction);
				double cosSqPhi = 1 - ((n * n) * (1 -(cosTheta * cosTheta))) / (nT * nT);
				double schlickReflectivity = 0.0;
				double schlickReflectivityTheta = 1.0;
				
				Color reflectColor = new Color();
				
		    	if (cosSqPhi > 0) // refraction
				{
		    		double cosPhi = Math.sqrt(cosSqPhi);
					schlickReflectivity = (nT - n) / (nT + n);
					schlickReflectivity *= schlickReflectivity;
					schlickReflectivityTheta = schlickReflectivity +  (1 - schlickReflectivity) * Math.pow(1 - cosTheta, 5);
					
					Vector3 refractedDir = new Vector3(iRec.frame.w);
					refractedDir.scale(cosTheta);
					refractedDir.add(ray.direction);
					refractedDir.scale(n / nT);
					refractedDir.scaleAdd(-cosPhi, iRec.frame.w);
					refractedRay.set(iRec.frame.o, refractedDir);
					refractedRay.makeOffsetRay();
					
					Ray outRay = new Ray();
					traceRefractedRay(scene, nT, refractedRay, 0, outRay);
					rayRadiance(scene, outRay, sampler, 0, reflectColor);
					outColor.scaleAdd(Math.max(0.0, 1 - schlickReflectivityTheta), reflectColor);
				}
				
		    	reflectColor.set(0,0,0);
				rayRadiance(scene, reflectedRay, sampler, 0, reflectColor);
		    	outColor.scaleAdd(Math.max(0.0, schlickReflectivityTheta), reflectColor);
			}
			
			//photon color
			Vector3 extents = scene.getBoundingBoxExtents();
			double sceneRadius = (extents.x > extents.y) ? ((extents.x > extents.z) ? extents.x : extents.z) 
									: ((extents.y > extents.z) ? extents.y : extents.z);
			//System.out.println( "Scene radius : " + sceneRadius);
			PhotonBoundingVolume.PhotonDistanceComparator comp = new PhotonBoundingVolume.PhotonDistanceComparator();
			PriorityQueue<PhotonBoundingVolume.HeapNode> globalHeap = new PriorityQueue<PhotonBoundingVolume.HeapNode>(PhotonMap.KNN, comp);
			PriorityQueue<PhotonBoundingVolume.HeapNode> causticHeap = new PriorityQueue<PhotonBoundingVolume.HeapNode>(PhotonMap.KNN, comp);
			this.photonMap.getFluxFromKNearest(iRec.frame.o, sceneRadius, globalHeap, causticHeap);
			Iterator<PhotonBoundingVolume.HeapNode> gIter = globalHeap.iterator();
			Iterator<PhotonBoundingVolume.HeapNode> cIter = causticHeap.iterator();
			//if(heap.size() < 2)
				//System.out.println(heap.size());
			
			//!!reverse incident directions for brdf?
			ray.direction.scale(-1.0);
			while (gIter.hasNext())
			{
				PhotonBoundingVolume.HeapNode heapNode = gIter.next();
				Color photonContribution = new Color(heapNode.photon.color);
				Color brdf = new Color();
				
				Vector3 incDir = new Vector3(heapNode.photon.incidentDirection);
				incDir.scale(-1.0);
				if(incDir.dot(iRec.frame.w) > 0)
				{
					iRec.surface.getMaterial().getBRDF(iRec).evaluate(iRec.frame, incDir, ray.direction, brdf);
//					if (heapNode.photon.isDirect)
//					{
						Color specBRDF = new Color();
						iRec.surface.getMaterial().getSpecular(iRec).evaluate(iRec.frame, incDir, ray.direction, specBRDF);
						brdf.add(specBRDF);
//					}
					photonContribution.scale(brdf);
					//TODO: Each radius or max radius
					photonContribution.scale(1 / (heapNode.distanceFromPoint * heapNode.distanceFromPoint * Math.PI));
					outColor.add(photonContribution);
				}
			}

			//ray.direction.scale(-1.0);
			
			while (cIter.hasNext())
			{
				PhotonBoundingVolume.HeapNode heapNode = cIter.next();
				Color photonContribution = new Color(heapNode.photon.color);
				Color brdf = new Color();
				Vector3 incDir = new Vector3(heapNode.photon.incidentDirection);
				incDir.scale(-1.0);
				if(incDir.dot(iRec.frame.w) > 0)
				{
					iRec.surface.getMaterial().getBRDF(iRec).evaluate(iRec.frame, incDir, ray.direction, brdf);
//					if (heapNode.photon.isDirect)
//					{
						Color specBRDF = new Color();
						iRec.surface.getMaterial().getSpecular(iRec).evaluate(iRec.frame, incDir, ray.direction, specBRDF);
						brdf.add(specBRDF);
//					}
					photonContribution.scale(brdf);
					//TODO: Each radius or max radius
					if (heapNode.distanceFromPoint > 0)
						photonContribution.scale(1 / (heapNode.distanceFromPoint * heapNode.distanceFromPoint * Math.PI));
					outColor.add(photonContribution);
				}
			}
			Color emitted = new Color();
			emittedRadiance(iRec, ray.direction, emitted);
			outColor.add(emitted);
			return;
			
		}			
		scene.getBackground().evaluate(ray.direction, outColor);
		*/
		
			}

		
	
	private void traceRefractedRay(Scene scene, double n, Ray inRay, int numInternalRefl, Ray outRay)
	{
		//Naming conventions from Shirley
		if (numInternalRefl > 6)
		{
			outRay = null;
			return;
		}
			
		double nT = 1.0;
		IntersectionRecord iRec = new IntersectionRecord();
		if (scene.getFirstIntersection(iRec, inRay))
		{
			iRec.frame.w.scale(-1.0);
			
			Vector3 reflectedDir = new Vector3(iRec.frame.w);
			reflectedDir.scale(-2.0 * (inRay.direction).dot(iRec.frame.w));
			reflectedDir.add(inRay.direction);
			
			Ray refractedRay = new Ray(inRay);
			refractedRay.direction.scale(-1.0);
			double cosTheta = iRec.frame.w.dot(refractedRay.direction);
			double cosSqPhi = 1 - ((n * n) * (1 -(cosTheta * cosTheta))) / (nT * nT);
		
			/*
			double schlickReflectivity = (nT - n) / (nT + n);
			schlickReflectivity *= schlickReflectivity;
			double schlickReflectivityTheta = schlickReflectivity +  (1 - schlickReflectivity) * Math.pow((1 - cosPhi), 5);*/
	        
	        if (cosSqPhi > 0) // refraction
			{
	        	double cosPhi = Math.sqrt(cosSqPhi);
				Vector3 refractedDir = new Vector3(iRec.frame.w);
				refractedDir.scale(cosTheta);
				refractedDir.add(inRay.direction);
				refractedDir.scale(n / nT);
				refractedDir.scaleAdd(-cosPhi, iRec.frame.w);
				outRay.set(iRec.frame.o, refractedDir);
				outRay.makeOffsetRay();
				return;
			}
			else
			{	
				//total internal reflection
				Ray reflectedRay = new Ray(iRec.frame.o, reflectedDir);
		        reflectedRay.makeOffsetRay();
				traceRefractedRay(scene, n, reflectedRay, numInternalRefl + 1, outRay);
			 }
		}		
	}
	
	protected void emittedRadiance(IntersectionRecord iRec, Vector3 dir, Color outColor) 
	{
		if(iRec.surface.getMaterial().isEmitter())
		{
			LuminaireSamplingRecord lRec = new LuminaireSamplingRecord();
			lRec.set(iRec);
            lRec.emitDir.set(dir);
            lRec.emitDir.scale(-1);
            iRec.surface.getMaterial().emittedRadiance(lRec, outColor);
			return;
		}
		outColor.set(0.0,0.0,0.0);
		
		// W4160 TODO (A)
		// If material is emitting, query it for emission in the relevant direction.
		// If not, the emission is zero.
		// This function should be called in the rayRadiance(...) method above
	}
	
	

}

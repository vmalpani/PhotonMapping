/*
 * Created on Nov 10, 2005
 * Copyright 2005 Program of Computer Grpahics, Cornell University
 */
package photon.accel;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.PriorityQueue;

import com.sun.org.apache.xerces.internal.impl.xpath.XPath.Axis;

import photon.Photon;
import photon.PhotonMap;

import ray.math.Point3;
import ray.misc.IntersectionRecord;
import ray.misc.Ray;
import ray.accel.AccelerationStructure;


/**
 * @author arbree
 * Nov 10, 2005
 * PhotonBoundingVolume.java
 * Copyright 2005 Program of Computer Graphics, Cornell University
 */
public class PhotonBoundingVolume implements AccelerationStructure {
	
	/** The maximum number of surfaces in a leaf node */
	public static final int MAX_PHOTONS_PER_LEAF = 20;
	
	/** The bounding box of this volume */
	protected final PhotonAxisAlignedBoundingBox box = new PhotonAxisAlignedBoundingBox();
	
	/** The surfaces contained in this bounding volume */
	protected ArrayList<Photon> photons = new ArrayList<Photon>();
	
	/** The children bounding volumes of this node */
	protected PhotonBoundingVolume left = null;
	protected PhotonBoundingVolume right = null;
	
	/** The depth of this node */
	protected final int depth;
	
	protected double splittingCoordinate;
	protected int splittingAxis = PhotonAxisAlignedBoundingBox.X;
	
	/**
	 * Private constructor used by split()
	 */
	private PhotonBoundingVolume(int inDepth) {
		depth = inDepth;
	}
	
	public PhotonAxisAlignedBoundingBox getBoundingBox() {
		return box;
	}
	
	/** Construct a bounding volume for the given surfaces, subdividing if necessary.
	 * @param inSurfaces
	 */
	public PhotonBoundingVolume(ArrayList<Photon> inPhotons) {
		
		depth = 0;
		
		//Add all the input surfaces to ourselves
		photons.addAll(inPhotons);
		
		growToHold();
		
		System.out.println("Volume contains "+inPhotons.size()+" objects.");
		
		//Divide if necessary
		split();
		
		System.out.println(box);
		
	}
	
	/**
	 * Grow the bounding volume to hold all the objects it encloses.
	 */
	private void growToHold() {
		
		//Grow our bounding box
		for (Iterator<Photon> iter = photons.iterator(); iter.hasNext();) {
			Photon currPhoton = (Photon) iter.next();
			currPhoton.addToBoundingBox(box);
		}
	}
	
	/**
	 * Split this bounding volume into two children
	 */
	private void split() {
		
		// If we are small enough, stop
		if(photons.size() < MAX_PHOTONS_PER_LEAF) {
			//vector.resize()
			photons.trimToSize();
			return;
		}
		
		//Create children
		left = new PhotonBoundingVolume(depth + 1);
		right = new PhotonBoundingVolume(depth + 1);
		
		//Break box along longest axis
		int axis = box.longestAxis();//depth % 3;
		Comparator<Photon> compare = null;
		switch(axis) {
		case PhotonAxisAlignedBoundingBox.X:
			compare = Photon.X_COMPARE;
			break;
		case PhotonAxisAlignedBoundingBox.Y:
			compare = Photon.Y_COMPARE;
			break;
		case PhotonAxisAlignedBoundingBox.Z:
			compare = Photon.Z_COMPARE;
			break;
		}
		
		//Sort the surfaces
		Collections.sort(photons, compare);
		
		//Put each half in the children
		
		ArrayList<Photon> leftList = new ArrayList<Photon>();
		
		List<Photon> firstHalf = photons.subList(0, photons.size()/2);
		leftList.addAll(firstHalf);
		firstHalf.clear();
		
		switch(axis) {
		case PhotonAxisAlignedBoundingBox.X:
			this.splittingCoordinate = ((Photon)(leftList.get(leftList.size() - 1))).getPosition().x
											+ ((Photon)(this.photons.get(0))).getPosition().x;
			this.splittingCoordinate /= 2.0;
			break;
		case PhotonAxisAlignedBoundingBox.Y:
			this.splittingCoordinate = ((Photon)(leftList.get(leftList.size() - 1))).getPosition().y
											+ ((Photon)(this.photons.get(0))).getPosition().y;
			this.splittingCoordinate /= 2.0;
			break;
		case PhotonAxisAlignedBoundingBox.Z:
			this.splittingCoordinate = ((Photon)(leftList.get(leftList.size() - 1))).getPosition().z
											+ ((Photon)(this.photons.get(0))).getPosition().z;
			this.splittingCoordinate /= 2.0;
			break;
		}
		
		//Set the object lists and clear ours
		
		left.photons = leftList;
		
		right.photons = this.photons;
		this.photons = null;
		
		this.splittingAxis = axis;
		
		//Grow children to fit
		left.growToHold();
		right.growToHold();
		
		left.split();
		right.split();
		
	}
	
	/**
	 * Set outRecord to the first intersection of ray with this bounding volume. Return true
	 * if there was an intersection and false otherwise. If no intersection was
	 * found outRecord is unchanged.
	 *
	 * @param outRecord the output IntersectionRecord
	 * @param ray the ray to intesect
	 * @return true if and intersection is found.
	 */
	
	// TODO: modify this function so as to get the k nearest neighbours
	
	public boolean getFirstIntersection(IntersectionRecord outRecord, Ray ray) {
		/*
		//Check that the ray intersects the box
		if(!box.intersect(ray))
			return false;
		
		//If we are a leaf, intersect our objects
		if(left == null && right == null) {
			
			//Find the first intersect by testing all surfaces
			double bestT = Double.MAX_VALUE;
			IntersectionRecord workRec = new IntersectionRecord();
			for (Iterator<Photon> iter = photons.iterator(); iter.hasNext();) {
				Photon currPhoton = (Photon) iter.next();
				if(currPhoton.intersect(workRec, ray) && workRec.t < bestT) {
					outRecord.set(workRec);
					bestT = workRec.t;
				
				}
			}
			
			return bestT != Double.MAX_VALUE;
			
		}
		
		//Check the left child
		IntersectionRecord leftRecord = new IntersectionRecord();
		if(left != null && left.getFirstIntersection(leftRecord, ray)) {
			
			//Shorten ray to hit point
			ray.end = leftRecord.t;
			
			//Check right child, the intersection must be closer than left
			IntersectionRecord rightRecord = new IntersectionRecord();
			if(right != null && right.getFirstIntersection(rightRecord, ray))
				outRecord.set(rightRecord);
			else outRecord.set(leftRecord);
			
			return true;
			
		}
		
		//Otherwise return the result of the right child
		return right != null && right.getFirstIntersection(outRecord, ray);
		*/
		return false;
	}
	
	public boolean getAnyIntersection(IntersectionRecord outRecord, Ray ray) {
		//could do this more efficiently.
		return getFirstIntersection(outRecord, ray);

	}
	
	public void getFluxFromKNearestRecursive(Point3 point, double queryRadiusSquare, PriorityQueue<HeapNode> heap)
	{
		if (this.right != null)
		{
			double distFromPlane = 0.0;
			switch(splittingAxis) {
			case PhotonAxisAlignedBoundingBox.X:
				distFromPlane = point.x - splittingCoordinate;
				break;
			case PhotonAxisAlignedBoundingBox.Y:
				distFromPlane = point.y - splittingCoordinate;
				break;
			case PhotonAxisAlignedBoundingBox.Z:
				distFromPlane = point.z - splittingCoordinate;
				break;
			}			
			
			if (distFromPlane < 0.0)
			{
				left.getFluxFromKNearestRecursive(point, queryRadiusSquare, heap);
				double distSquare = distFromPlane * distFromPlane;

				// (Double.compare(distSquare, queryRadiusSquare) > 0)
				if (distSquare < queryRadiusSquare)
				{
					/*System.out.println("distSquare : " + distSquare);
					System.out.println("queryRadius : " + queryRadiusSquare);*/

					right.getFluxFromKNearestRecursive(point, queryRadiusSquare, heap);
				}
			}
			else
			{
				double distSquare = distFromPlane * distFromPlane;
				right.getFluxFromKNearestRecursive(point, queryRadiusSquare, heap);
				//if (Double.compare(distSquare, queryRadiusSquare) < 0)
				if (distSquare < queryRadiusSquare)
				{
					/*System.out.println("distSquare : " + distSquare);
					System.out.println("queryRadius : " + queryRadiusSquare);*/

					left.getFluxFromKNearestRecursive(point, queryRadiusSquare, heap);
				}
			}
		}
		else if (null != this.photons)
		{
			//nbPhotonsVisited[0] = nbPhotonsVisited[0] + 1;//this.photons.size();

			for (int i = 0; i < this.photons.size(); i++)
			{
				double distFromPoint = point.distanceSquared(this.photons.get(i).getPosition());
				if (Double.compare(distFromPoint, queryRadiusSquare) < 0)
				{
					HeapNode h = new HeapNode(distFromPoint, this.photons.get(i));
					if (heap.size() < PhotonMap.KNN)
						heap.add(h);
					else if (heap.peek().distanceFromPoint > distFromPoint)
					{
						heap.remove();
						heap.add(h);
					}
				}
			}
		}	
	}
	
	public static class HeapNode
	{
		public double distanceFromPoint = 0.0;
		public Photon photon = null;
		
		HeapNode(double d, Photon p)
		{
			this.distanceFromPoint = d;
			this.photon = p;
		}
	}
	
	public static class PhotonDistanceComparator implements Comparator<PhotonBoundingVolume.HeapNode>
	{
		public int compare(PhotonBoundingVolume.HeapNode n0, PhotonBoundingVolume.HeapNode n1) {
			// TODO Auto-generated method stub
			return Double.compare(n1.distanceFromPoint, n0.distanceFromPoint);
		}
	}
	
}

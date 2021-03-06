package ray.surface;

import ray.accel.AxisAlignedBoundingBox;
import ray.material.Material;
import ray.math.Frame3;
import ray.math.Geometry;
import ray.math.Point2;
import ray.math.Point3;
import ray.math.Vector3;
import ray.misc.IntersectionRecord;
import ray.misc.LuminaireSamplingRecord;
import ray.misc.Ray;

/**
 * Represents a sphere as a center and a radius.
 *
 * @author ags
 */
public class Sphere extends Surface {
    
    /** Material for this sphere. */
    protected Material material = Material.DEFAULT_MATERIAL;
    
    /** The center of the sphere. */
    protected final Point3 center = new Point3();
    
    /** The radius of the sphere. */
    protected double radius = 1.0;
    
    private AxisAlignedBoundingBox boundingBox = new AxisAlignedBoundingBox();
    
    /**
     * Default constructor, creates a sphere at the origin with radius 1.0
     */
    public Sphere() { 
    	boundingBox = new AxisAlignedBoundingBox();
        this.boundingBox.add(center.x - radius, center.y - radius, center.z - radius);
        this.boundingBox.add(center.x + radius, center.y + radius, center.z + radius);
    }
    
    /**
     * The explicit constructor. This is the only constructor with any real code
     * in it. Values should be set here, and any variables that need to be
     * calculated should be done here.
     *
     * @param newCenter The center of the new sphere.
     * @param newRadius The radius of the new sphere.
     * @param newMaterial The material of the new sphere.
     */
    public Sphere(Vector3 newCenter, double newRadius, Material newMaterial) {        
        material = newMaterial;
        center.set(newCenter);
        radius = newRadius;
        updateArea();
        this.boundingBox.add(center.x - radius, center.y - radius, center.z - radius);
        this.boundingBox.add(center.x + radius, center.y + radius, center.z + radius);
    }
    
    public void updateArea() {
    	area = 4 * Math.PI * radius*radius;
    	oneOverArea = 1. / area;
    }
    
    /**
     * @see ray1.surface.Surface#getMaterial()
     */
    public Material getMaterial() {
        return material;
    }
    
    /**
     * @see ray1.surface.Surface#setMaterial(ray1.material.Material)
     */
    public void setMaterial(Material material) {
        this.material = material;
    }
    
    /**
     * Returns the center of the sphere in the input Point3
     * @param outPoint output space
     */
    public void getCenter(Point3 outPoint) {        
        outPoint.set(center);    
        
        
    }
    
    /**
     * @param center The center to set.
     */
    public void setCenter(Point3 center) {        
        this.center.set(center);  
        boundingBox = new AxisAlignedBoundingBox();
        this.boundingBox.add(center.x - radius, center.y - radius, center.z - radius);
        this.boundingBox.add(center.x + radius, center.y + radius, center.z + radius);
    }
    
    /**
     * @return Returns the radius.
     */
    public double getRadius() {
        return radius;
    }
    
    /**
     * @param radius The radius to set.
     */
    public void setRadius(double radius) {
        this.radius = radius;
        boundingBox = new AxisAlignedBoundingBox();
        this.boundingBox.add(center.x - radius, center.y - radius, center.z - radius);
        this.boundingBox.add(center.x + radius, center.y + radius, center.z + radius);
        updateArea();
    }
    
    public boolean chooseSamplePoint(Point3 p, Point2 seed, LuminaireSamplingRecord lRec) {
        Geometry.squareToSphere(seed, lRec.frame.w);
        lRec.frame.o.set(center);
        lRec.frame.o.scaleAdd(radius, lRec.frame.w);
        lRec.frame.initFromW();
        lRec.pdf = oneOverArea;
        lRec.emitDir.sub(p, lRec.frame.o);
        return true;
    }
        
    /**
     * @see ray1.surface.Surface#intersect(ray1.misc.IntersectionRecord,
     *      ray1.misc.Ray)
     */
    public boolean intersect(IntersectionRecord outRecord, Ray ray) {
        // W4160 TODO (A)
    	// In this function, you need to test if the given ray intersect with a sphere.
    	// You should look at the intersect method in other classes such as ray.surface.Triangle.java
    	// to see what fields of IntersectionRecord should be set correctly.
    	// The following fields should be set in this function
    	//     IntersectionRecord.t
    	//     IntersectionRecord.frame
    	//     IntersectionRecord.surface
    	//
    	// Note: Although a ray is conceptually a infinite line, in practice, it often has a length,
    	//       and certain rendering algorithm relies on the length. Therefore, here a ray is a 
    	//       segment rather than a infinite line. You need to test if the segment is intersect
    	//       with the sphere. Look at ray.misc.Ray.java to see the information provided by a ray.
    	Vector3 ce = new Vector3();
    	ce.sub(ray.origin, this.center);
    	double dDotd = ray.direction.squaredLength();
    	double dDotce = ray.direction.dot(ce); 
    	double discriminant = dDotce * dDotce - dDotd * (ce.squaredLength() - this.radius * this.radius);
    	if (discriminant < 0)
    	{
    		//No intersection
    		return false;
    	}
    	double t = (-dDotce - Math.sqrt(discriminant)) / dDotd;
    	if (t < ray.start || t > ray.end)
    	{
    		t = (-dDotce + Math.sqrt(discriminant)) / dDotd;
    		if (t < ray.start || t > ray.end)
    			return false;
    		Point3 intersectionPoint = new Point3(ray.origin);
	    	intersectionPoint.scaleAdd(t, ray.direction);
	    	Vector3 normal = new Vector3();
	    	normal.sub(intersectionPoint, this.center);
	    	
	    	outRecord.frame.set(intersectionPoint, new Vector3(1, 0, 0), new Vector3(0, 1, 0), normal);
	    	outRecord.frame.initFromW();
	    	
	    	outRecord.t = t;
	    	outRecord.surface = this;
    	}
    	
    	else
    	{
	    	Point3 intersectionPoint = new Point3(ray.origin);
	    	intersectionPoint.scaleAdd(t, ray.direction);
	    	Vector3 normal = new Vector3();
	    	normal.sub(intersectionPoint, this.center);
	    	
	    	outRecord.frame.set(intersectionPoint, new Vector3(1, 0, 0), new Vector3(0, 1, 0), normal);
	    	outRecord.frame.initFromW();
	    	
	    	outRecord.t = t;
	    	outRecord.surface = this;
    	}
        
    	return true;
    }
    
    /**
     * @see Object#toString()
     */
    public String toString() {
        return "sphere " + center + " " + radius + " " + material + " end";
    }
    
    /**
     * @see ray1.surface.Surface#addToBoundingBox(ray1.surface.accel.AxisAlignedBoundingBox)
     */
    public void addToBoundingBox(AxisAlignedBoundingBox inBox) {
        inBox.add(center.x - radius, center.y - radius, center.z - radius);
        inBox.add(center.x + radius, center.y + radius, center.z + radius);        
    }

    public AxisAlignedBoundingBox getBoundingBox() {
    	return this.boundingBox;
	}
    
}

package photon;

import java.util.Comparator;

import photon.accel.PhotonAxisAlignedBoundingBox;
import ray.math.Point3;
import ray.math.Vector3;
import ray.misc.Color;

public class Photon {
	
	public static final Comparator<Photon> X_COMPARE = new PhotonComparator(PhotonComparator.X_AXIS);
	public static final Comparator<Photon> Y_COMPARE = new PhotonComparator(PhotonComparator.Y_AXIS);
	public static final Comparator<Photon> Z_COMPARE = new PhotonComparator(PhotonComparator.Z_AXIS);
	public Point3 position = new Point3();
	public Color color = new Color();
	public Vector3 incidentDirection = new Vector3();
	public boolean isDirect = false;
	
	Photon (){
		
	}
	
	public void addToBoundingBox(PhotonAxisAlignedBoundingBox inBox)
	{
		inBox.add(this.position);
	}

	private static class PhotonComparator implements Comparator<Photon> {
        
        //Constants defining the axis directions
        public static final int X_AXIS = 0;
        public static final int Y_AXIS = 1;
        public static final int Z_AXIS = 2;
        
        //The axis of this comparator
        protected int axis;
        
        /*//Stores the positions of the objects currently being compared
        private Point3 c1 = new Point3();
        private Point3 c2 = new Point3();*/
        
        /**
         * Constructs a comparator for a certain axis
         * @param inAxis
         */
        public PhotonComparator(int inAxis) {
            
            axis = inAxis;
            
        }
        
        /**
         * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
         */
        public int compare(Photon s1, Photon s2) {
            
            Point3 c1 = s1.getPosition();
            Point3 c2 = s2.getPosition();
            
            //Get the direction to compare
            double d1 = 0;
            double d2 = 0;
            switch (axis) {
            case X_AXIS:
                d1 = c1.x;
                d2 = c2.x;
                break;
            case Y_AXIS:
                d1 = c1.y;
                d2 = c2.y;
                break;
            case Z_AXIS:
                d1 = c1.z;
                d2 = c2.z;
                break;
            }
            
            return Double.compare(d1, d2);
            
        }
    }

	public Point3 getPosition() {
		return this.position;
	}
}

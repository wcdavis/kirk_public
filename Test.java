import com.vividsolutions.jts.geom.Coordinate;

import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.Polygon;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.triangulate.DelaunayTriangulationBuilder;
import com.vividsolutions.jts.*;

import org.jgrapht.experimental.dag.DirectedAcyclicGraph;
import org.jgrapht.graph.SimpleDirectedGraph;
import org.jgrapht.graph.SimpleGraph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.alg.CycleDetector;

import java.util.List;
import java.util.Collection;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.HashSet;
public class Test
{
    private Geometry finalTriangle;
    private SimpleDirectedGraph<Geometry, DefaultEdge> dag;
    private SimpleGraph<Coordinate, DefaultEdge> graph;
    private Set<Geometry> validTriangles;
    private CycleDetector<Geometry, DefaultEdge> cd;

    public Test(Coordinate[] points)
    {
		/* Initialize local variables */
		int i;
		int pointCounter  = 0;
		int regionCounter = 0;
		DelaunayTriangulationBuilder dt = new DelaunayTriangulationBuilder();
        GeometryFactory fact = new GeometryFactory();
        HashMap<Coordinate, List<Geometry>> coordinateToGeometry = new HashMap<Coordinate, List<Geometry>>();
        HashMap<Coordinate, List<Coordinate>> coordinateToNeighbors = new HashMap<Coordinate, List<Coordinate>>();

		/* Initialize private variables */
		dag            = new SimpleDirectedGraph<Geometry, DefaultEdge>(DefaultEdge.class);
		graph          = new SimpleGraph<Coordinate, DefaultEdge>(DefaultEdge.class);
		validTriangles = new HashSet<Geometry>();

        cd = new CycleDetector<>(dag);
		/* Create outer bounding triangle to contain all the regions */
		Coordinate p1       = new Coordinate(-1000000, -1000000); 
		Coordinate p2       = new Coordinate(0.0,                  1000000); 
		Coordinate p3       = new Coordinate(1000000,  -1000000); 

        Coordinate[] boundaryPointsArray = new Coordinate[4];
        boundaryPointsArray[0] = p1;
        boundaryPointsArray[1] = p2;
        boundaryPointsArray[2] = p3;
        boundaryPointsArray[3] = p1;

        ArrayList<Coordinate> boundaryPoints = new ArrayList<Coordinate>();
        boundaryPoints.add(p1);
        boundaryPoints.add(p2);
        boundaryPoints.add(p3);

        Polygon finalTrianglePoly = fact.createPolygon(boundaryPointsArray);
        finalTriangle = finalTrianglePoly.convexHull();

        Set<Coordinate> allPoints = new HashSet<Coordinate>();
        for(Coordinate p : points)
        {
            allPoints.add(p);
        }

        /* We get the triangulation of just the points so that we can log
         * which triangles are from the original shape */
		dt.setSites(allPoints);
        Geometry g = dt.getTriangles(fact);
        List<Geometry> triangles = new ArrayList<Geometry>();
        for (i = 0; i < g.getNumGeometries(); i++)
        {
            triangles.add(g.getGeometryN(i));
        }

        /* add triangles to the valid list */
        for (Geometry t : triangles)
		{
            validTriangles.add(t);
		}

        /* Add boundary points and triangulate between convex hull and shell */
        dt = new DelaunayTriangulationBuilder();
        allPoints.addAll(boundaryPoints);
        dt.setSites(allPoints);

        g = dt.getTriangles(fact);
        triangles = new ArrayList<Geometry>();
        for (i = 0; i < g.getNumGeometries(); i ++)
        {
            triangles.add(g.getGeometryN(i));
        }

        /* Create graph that we later use to find independent sets */
        for (Geometry t : triangles)
		{
            Coordinate[] coordinates = t.getCoordinates();
            graph.addVertex(coordinates[0]);
            graph.addVertex(coordinates[1]);
            graph.addVertex(coordinates[2]);
			graph.addEdge(coordinates[0], coordinates[1]);
			graph.addEdge(coordinates[0], coordinates[2]);
			graph.addEdge(coordinates[1], coordinates[2]);

            /* add triangles to the digraph as vertices for now */
            dag.addVertex(t);

            for (i = 0; i < 3; i++)
            {
                List<Geometry> geoList;
                List<Coordinate> neighborList; 


                if (coordinateToGeometry.get(coordinates[i]) != null)
                    geoList = coordinateToGeometry.get(coordinates[i]);

                else 
                    geoList = new ArrayList<Geometry>();

                geoList.add(t);
                coordinateToGeometry.put(coordinates[i], geoList);

                if (coordinateToNeighbors.get(coordinates[i]) != null)
                    neighborList = coordinateToNeighbors.get(coordinates[i]);

                else 
                    neighborList = new ArrayList<Coordinate>();

                for (int j = 0; j < 4; j++)
                {
                    Coordinate c = coordinates[j];
                    if(c.equals(coordinates[i]))
                        continue;
                    neighborList.add(c);
                }

                coordinateToNeighbors.put(coordinates[i], neighborList);
            }
        }

        do
        {

            /* Create independent set. Use two lists to avoid modifying list we're iterating over */
			ArrayList<Coordinate> independentSetCandidates = new ArrayList<Coordinate>();
			ArrayList<Coordinate> independentSet            = new ArrayList<Coordinate>();

            Set<Coordinate> vertices = graph.vertexSet();
			for (Coordinate p : vertices)
			{

				if (graph.degreeOf(p) <= 8 && !boundaryPoints.contains(p))
                {

					independentSetCandidates.add(p);
                }
			}

			for (Coordinate u : independentSetCandidates)
			{
                boolean add = true;
                for (Coordinate v : independentSet)
                {
                    if (graph.containsEdge(u,v) || graph.containsEdge(v,u))
                        add = false;

                    if (v.equals(u))
                        add = false;
                }

                if (add)
                    independentSet.add(u);

			}

            

            for(Coordinate p : independentSet)
            {
            /* update neighbor and affected triangles list */
            dt = new DelaunayTriangulationBuilder();
            dt.setSites(graph.vertexSet());
            triangles = new ArrayList<Geometry>();
            g = dt.getTriangles(fact);
            for (i = 0; i < g.getNumGeometries(); i ++)
            {
                triangles.add(g.getGeometryN(i));
            }

            coordinateToNeighbors = new HashMap<Coordinate, List<Coordinate>>(); 
            coordinateToGeometry  = new HashMap<Coordinate, List<Geometry>>(); 

            for (Geometry t : triangles)
            {
                List<Geometry> geoList;
                List<Coordinate> neighborList; 
                Coordinate[] coordinates = t.getCoordinates();


                for (i = 0; i < 3; i++)
                {
                        if (coordinateToGeometry.get(coordinates[i]) != null)
                            geoList = coordinateToGeometry.get(coordinates[i]);

                        else 
                            geoList = new ArrayList<Geometry>();

                        geoList.add(t);
                        coordinateToGeometry.put(coordinates[i], geoList);

                        if (coordinateToNeighbors.get(coordinates[i]) != null)
                            neighborList = coordinateToNeighbors.get(coordinates[i]);

                        else 
                            neighborList = new ArrayList<Coordinate>();

                        for (int j = 0; j < 4; j++)
                        {
                            Coordinate c = coordinates[j];
                            if(c.equals(coordinates[i]))
                                continue;
                            neighborList.add(c);
                        }

                        coordinateToNeighbors.put(coordinates[i], neighborList);
                    }
            }
                ArrayList<Coordinate> pCollection = new ArrayList<Coordinate>();
                pCollection.add(p);

                Collection<DefaultEdge> edges = graph.edgesOf(p);
                List<Coordinate> neighbors = new ArrayList<Coordinate>();
                List<Geometry> oldTriangles = new ArrayList<Geometry>();
                List<Geometry> newTriangles = new ArrayList<Geometry>();

                neighbors = coordinateToNeighbors.get(p);
                dt = new DelaunayTriangulationBuilder();
                dt.setSites(neighbors);
                g = dt.getTriangles(fact);
                for (i = 0; i < g.getNumGeometries(); i ++)
                {
                    newTriangles.add(g.getGeometryN(i));
                }

                oldTriangles = coordinateToGeometry.get(p);
                oldTriangles.removeAll(newTriangles);

                
                if (newTriangles.size() == 1)
                    finalTriangle = newTriangles.get(0);

                /* Update undirected graph */
                graph.removeVertex(p);

                /* Update DAG */
                for (Geometry newT : newTriangles)
                {
                    dag.addVertex(newT);

                    for (Geometry oldT : oldTriangles)
                    {
                        dag.addVertex(oldT);
                        dag.addEdge(newT, oldT);
                    }
                }
            }

        } while(graph.vertexSet().size() > 3);
    }

    public Geometry lookup(Coordinate p)
    {
        GeometryFactory fact = new GeometryFactory();
        Point point = fact.createPoint(p);

		Geometry curr = finalTriangle;

        if (cd.detectCycles())

		if (!curr.contains(point))
		{
			return null;
		}

		while(!dag.outgoingEdgesOf(curr).isEmpty())
		{
			for(DefaultEdge e : dag.outgoingEdgesOf(curr))
			{
				Geometry t = dag.getEdgeTarget(e);
				if (t.contains(point))
				{
					curr = t;
					break;
				}
			}
		}

        for (Geometry g : validTriangles)
        {
            if (g.contains(point))
                System.out.println("************ " + g + " *************");
            if (g.equals(curr))
            {
                assert(g.contains(point));
                return curr;
            }
        }
		return curr;

    }

    public static void main (String[] args)
    {
        
        int numPoints = Integer.parseInt(args[0]);
        Coordinate[] points = new Coordinate[numPoints];
        Coordinate c1 = new Coordinate(0.0, 0.0);
        Coordinate c2 = new Coordinate(1.0, 0.0);
        Coordinate c3 = new Coordinate(1.0, 1.0);
        Coordinate c4 = new Coordinate(0.0, 1.0);
        points[0] = c1;
        points[1] = c2;
        points[2] = c3;
        points[3] = c4;
        

        for (int i = 4; i < numPoints; i++)
        {
            double x = Math.random();
            double y = Math.random();
            Coordinate p = new Coordinate(x,y);
            points[i] = p;
        }


        double start = System.currentTimeMillis();
        Test t = new Test(points);
        double end = System.currentTimeMillis();
        System.out.println(end - start);
        start = System.currentTimeMillis();
        t.lookup(new Coordinate(Math.random(),Math.random()));
        end = System.currentTimeMillis();
        System.out.println(end - start);
        start = System.currentTimeMillis();
        t.lookup(new Coordinate(0.1,0.5));
        end = System.currentTimeMillis();
        System.out.println(end - start);
    }
}


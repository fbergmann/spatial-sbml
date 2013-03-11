import java.awt.BorderLayout;
import java.awt.EventQueue;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.Timer;
import org.spatial.SpatialSimulator;

/**
 * Simple example of calling the spatialSBML Java API
 * 
 * @author Frank T. Bergmann <fbergman@caltech.edu>
 * 
 */
public class SpatialWindow
{

    // only need to change these constants
    /**
     * The step size used by the simulator
     */
    final double        stepSize        = 0.08;
    /**
     * The file name to be loaded and simulated
     */
    final String        fileName        = "D:/Development/spatial-sbml/examples/turing.xml";
    /**
     * SBML Id of species to observe
     */
    final String        observedSpecies = "species_2";
    /**
     * Maximum concentration of observed species
     */
    final static double max             = 1.5f;

    /**
     * Scale factor to scale pixels by
     */
    final int           scale           = 4;

    /**
     * Number of milliseconds the timer sleeps between simulation steps.
     */
    final int           speed           = 5;

    // should not need to modify anything below

    /**
     * Not quite sure how to render the image properly, so I
     * basically override a panel and keep painting.
     * 
     */
    private class ImagePanel
        extends JPanel
    {

        /**
         * The image to paint
         */
        private Image image;


        /**
         * Construct panel
         */
        public ImagePanel()
        {
            super();
        }


        /**
         * Set the image
         * 
         * @param image
         *            the image to set
         */
        public void setImage(BufferedImage image)
        {
            this.image = image;
            repaint();
        }


        /**
         * And draw the image ... if not null
         */
        @Override
        public void paintComponent(Graphics g)
        {
            super.paintComponent(g);
            if (image != null)
            {
                g.drawImage(image, 0, 0, this);
            }
        }

    }

    /**
     * Basic timer, that steps through the simulation.
     * 
     */
    private final class SimulationTimer
        implements ActionListener
    {
        double time  = 0;
        long   count = 0;


        /**
         * Constructor with filename
         */
        private SimulationTimer(String fileName)
        {
            sim = new SpatialSimulator();
            sim.initFromFile(fileName, 101, 101);
        }


        /**
         * perform simulation steps, and call update after 10 steps
         */
        @Override
        public void actionPerformed(ActionEvent e)
        {
            count++;
            time = sim.oneStep(time, stepSize);

            if (count % 10 == 0)

            {
                drawImage(sim, time);
            }

        }
    }

    private JFrame frmSpatialTestjava;


    /**
     * Launch the application.
     */
    public static void main(String[] args)
    {
        /**
         * Load the native library
         */
        System.loadLibrary("SpatialJava");

        EventQueue.invokeLater(new Runnable() {
            @Override
            public void run()
            {
                try
                {
                    SpatialWindow window = new SpatialWindow();
                    window.frmSpatialTestjava.setVisible(true);
                }
                catch (Exception e)
                {
                    e.printStackTrace();
                }
            }
        });
    }


    /**
     * Create the application.
     */
    public SpatialWindow()
    {
        initialize();
    }

    private SpatialSimulator sim;
    private Timer            timer;
    private ImagePanel       panel;


    /**
     * Initialize the contents of the frame.
     */
    private void initialize()
    {

        frmSpatialTestjava = new JFrame();
        frmSpatialTestjava.setTitle("Spatial Test (Java API)");
        frmSpatialTestjava.setBounds(100, 100, 518, 555);
        frmSpatialTestjava.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frmSpatialTestjava.getContentPane().setLayout(new BorderLayout(0, 0));

        panel = new ImagePanel();
        frmSpatialTestjava.getContentPane().add(panel);

        lblTime = new JLabel("");
        frmSpatialTestjava.getContentPane().add(lblTime, BorderLayout.SOUTH);

        JMenuBar menuBar = new JMenuBar();
        frmSpatialTestjava.setJMenuBar(menuBar);

        JMenu mnNewMenu = new JMenu("File");
        menuBar.add(mnNewMenu);

        JMenuItem mnOpen = new JMenuItem("Open");
        mnOpen.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e)
            {
                if (timer != null)
                {
                    timer.stop();
                }

                timer = new Timer(speed, new SimulationTimer(fileName));
                timer.setInitialDelay(speed);
                timer.start();

            }
        });
        mnNewMenu.add(mnOpen);

        JMenuItem mnExit = new JMenuItem("Exit");
        mnExit.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e)
            {
                System.exit(0);
            }
        });
        mnNewMenu.add(mnExit);
    }


    /**
     * Draw the image by getting the concentrations of the observed variable
     * from the
     * given simulator.
     * 
     * @param sim2
     *            the simulator
     * @param time
     *            current step
     */
    protected synchronized void drawImage(SpatialSimulator sim2, double time)

    {
        // long start = System.currentTimeMillis();
        int xDim = sim2.getXDim();
        int yDim = sim2.getYDim();
        final BufferedImage image = new BufferedImage(xDim * scale, yDim
            * scale, java.awt.image.BufferedImage.TYPE_INT_RGB);
        for (int y = 0; y < yDim; y += 1)
        {
            for (int x = 0; x < xDim; x += 1)
            {
                double current = sim.getVariableAt(observedSpecies, x, y);

                int rgb = getRgbForConcentration(current);
                if (rgb != 0)
                {
                    image.setRGB(x * scale, y * scale, rgb);
                }
            }
        }
        panel.setImage(image);
        lblTime.setText(String.format("Time: %f", time));
        // System.out.println("Drawing tool: " + (System.currentTimeMillis() -
        // start));

    }

    private JLabel lblTime;


    /**
     * Utility function mapping the current concentration to an RGB value.
     * 
     * @param current
     *            current concentration
     * @return rgb value for that concentration
     */
    private int getRgbForConcentration(double current)
    {
        double maxV = Math.max(current, max);
        int color = (int) ((current / maxV) * (256 - 1));
        int r = 0;// red component 0...255
        int g = color;// green component 0...255
        int b = 0;// blue component 0...255
        return (r << 16) | (g << 8) | b;

    }
}

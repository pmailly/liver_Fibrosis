package Heart_fibrosis_CD_tools;




import ij.ImagePlus;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.filter.ParticleAnalyzer;
import ij.process.ImageProcessor;
import java.awt.Rectangle;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

        
 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose dots_Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author phm
 */
public class Heart_fibrosist_tools {    

    /**
     * Find serie from roi filename
     * @param files
     * @param fileName
     * @return 
     */
    public static String find_roi(String[] files, String fileName) {
        String series = "";
        for (String f : files) {
            if (f.matches(fileName+"#.*.zip"))
                series = f.split("#")[1];
	}
	return (series);
    }
    
    /**
     * Find color matrix from roi
     * @param img
     * @param rgbOD
     * @param series
     * @param roi
     * @return 
     */
    public static double[] findColorMatrix(ImagePlus img, double[] rgbOD, Roi roi) {
	int  x, y, px, py, w, h, p;
        double log255 = Math.log(255.0);
        ImageProcessor ip = img.getProcessor();
        img.setRoi(roi);
        img.updateAndDraw();
        Rectangle bounds = roi.getBounds();
        px = bounds.x;
        py = bounds.y;
        h = bounds.height;
        w = bounds.width;
        
	for (x = px; x < (px + w); x++){
            for(y = py; y < (py + h); y++){
                p = ip.getPixel(x,y);
              // rescale to match original paper values
              rgbOD[0] = rgbOD[0]+ (-((255.0*Math.log(((double)((p & 0xff0000)>>16)+1)/255.0))/log255));
              rgbOD[1] = rgbOD[1]+ (-((255.0*Math.log(((double)((p & 0x00ff00)>> 8) +1)/255.0))/log255));
              rgbOD[2] = rgbOD[2]+ (-((255.0*Math.log(((double)((p & 0x0000ff))        +1)/255.0))/log255));
            }
        }
        rgbOD[0] = rgbOD[0] / (w*h);
        rgbOD[1] = rgbOD[1] / (w*h);
        rgbOD[2] = rgbOD[2] / (w*h);
        return(rgbOD);
    }
  
    public static double analyse(ImagePlus img, String path ) {
        ResultsTable rt = new ResultsTable();
        ParticleAnalyzer ana = new ParticleAnalyzer(ParticleAnalyzer.CLEAR_WORKSHEET+ParticleAnalyzer.SHOW_ROI_MASKS, Measurements.AREA+Measurements.LIMIT, rt, 20000, Double.POSITIVE_INFINITY);
        ana.analyze(img);
        double area = rt.getValue("Area", 0);
        System.out.println("red ="+area);
        img.close();
        ImagePlus imgRedCount = ana.getOutputImage();
        // Save
        FileSaver imgSave = new FileSaver(imgRedCount);
        imgSave.saveAsTiff(path);
        imgRedCount.close();
        return(area);
    }
    
    
     /**
     * 
     * @param outDirResults
     * @param resultsFileName
     * @param header
     * @return 
     * @throws java.io.IOException 
     */
    public static BufferedWriter writeHeaders(String outDirResults, String resultsFileName, String header) throws IOException {
        FileWriter FileResults = new FileWriter(outDirResults + resultsFileName, false);
        BufferedWriter outPutResults = new BufferedWriter(FileResults); 
        outPutResults.write(header);
        outPutResults.flush();
        return(outPutResults);
    }        
}

package Liver_fibrosis_CD_tools;




import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.io.Opener;
import ij.measure.Calibration;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.ParticleAnalyzer;
import ij.plugin.frame.RoiManager;
import ij.process.ImageConverter;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij.clearcl.ClearCLImage;
import net.haesleinhuepf.clij2.CLIJ2;

        
 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose dots_Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author phm
 */
public class Liver_fibrosis_tools { 
    
    public static double tileSize = 6000;
    public static boolean saveMask = false;
    private static String[] color_staining = {"Sirius red", "Tunel"};
    public static String color_meth = "";
    public static CLIJ2 clij2 = CLIJ2.getInstance();
    
    
    
     /**
     * check  installed modules
     * @return 
     */
    public static boolean checkInstalledModules() {
        // check install
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("net.haesleinhuepf.clij2.CLIJ2");
        } catch (ClassNotFoundException e) {
            IJ.log("CLIJ not installed, please install from update site");
            return false;
        }
        try {
            loader.loadClass("mcib3d.geom");
        } catch (ClassNotFoundException e) {
            IJ.log("3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        return true;
    }

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
     * Dialog 

     */
    public static String dialog() {
        String dir = "";
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.addDirectoryField("Choose Directory Containing Image Files : ", "");
        gd.addChoice("Staining method : ", color_staining, color_staining[0]);
        gd.addNumericField("Tile size : ", tileSize, 0);
        gd.addCheckbox(" Save mask", saveMask);
        gd.showDialog();
        dir = gd.getNextString();
        color_meth = gd.getNextChoice();
        tileSize = gd.getNextNumber();
        saveMask = gd.getNextBoolean();
        return(dir);
    }
    
     /**  
     * median 2D box filter
     * Using CLIJ2
     * @param imgCL
     * @param sizeX
     * @param sizeY
     * @param sizeZ
     * @return imgOut
     */ 
    public static ClearCLBuffer medianFilter(ClearCLBuffer imgCL, double sizeX, double sizeY) {
        ClearCLBuffer imgOut = clij2.create(imgCL);
        clij2.median2DBox(imgCL, imgOut, sizeX, sizeY);
        clij2.release(imgCL);
        return(imgOut);
    }
    
      
    
     /**  
     * exponentiel box filter
     * Using CLIJ2
     * @param imgCL
     * @return imgOut
     */ 
    public static ClearCLBuffer clij_filter(ClearCLBuffer imgCL, String filter) {
        ImagePlus img32 = clij2.pull(imgCL);
        new ImageConverter(img32).convertToGray32();
        ClearCLBuffer imgIn = clij2.push(img32);
        ClearCLBuffer imgOut = clij2.create(imgIn);
        if (filter.equals("exp"))
            clij2.exponential(imgIn, imgOut);
        else
            clij2.logarithm(imgIn, imgOut);    
        
        clij2.release(imgCL);
        clij2.release(imgIn);
        img32.close();
        return(imgOut);
    }
    
    
    
  /**
   * Open
   * USING CLIJ2
   * @param imgCL
   * @return imgCLOut
   */
    private static ClearCLBuffer open(ClearCLBuffer imgCL, int iter) {
        ClearCLBuffer imgCLOut = clij2.create(imgCL);
        clij2.openingBox(imgCL, imgCLOut, iter);
        clij2.release(imgCL);
        return(imgCLOut);
    }
    
    /**
     * Threshold 
     * USING CLIJ2
     * @param imgCL
     * @param thMed
     */
    public static ClearCLBuffer threshold(ClearCLBuffer imgCL, String thMed) {
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        clij2.release(imgCL);
        return(imgCLBin);
    }
    
    /**
     * 
     * @param img
     * @param imgFibrosis
     * @param th
     * @param filter
     * @param open
     * @param min
     * @param max
     * @param path
     * @return 
     */
    
     public static double clij_analyse(ImagePlus img, ImagePlus imgFibrosis, String th, String filter, int open, double min, double max, String path) {
        double fibArea = 0;
        double circMax = 1;
        ResultsTable rt = new ResultsTable();
        Calibration cal = imgFibrosis.getCalibration();
        
        ClearCLBuffer imgCL = clij2.push(imgFibrosis);
        imgFibrosis.close();
        ClearCLBuffer imgCLmed = medianFilter(imgCL, 4, 4);
        ClearCLBuffer imgCLfilter;
        if (filter.equals("exp"))
            imgCLfilter = clij_filter(imgCLmed, "exp");
        else {
            imgCLfilter = clij_filter(imgCLmed, "log");
            circMax = 0.7;
        }
        ClearCLBuffer imgCLbin = threshold(imgCLfilter, th);
        ImagePlus imgBin;
        if (open != 0) {
            ClearCLBuffer imgCLopen = open(imgCLbin, open);
            imgBin = clij2.pull(imgCLopen);
            clij2.release(imgCLopen);
        }
        else 
            imgBin = clij2.pull(imgCLbin);
        clij2.release(imgCLbin);
        IJ.run(imgBin, "8-bit", "");
        imgBin.setCalibration(cal);
        ParticleAnalyzer ana = new ParticleAnalyzer(ParticleAnalyzer.CLEAR_WORKSHEET+ParticleAnalyzer.SHOW_ROI_MASKS,
                Measurements.AREA+Measurements.MEAN, rt, min, max, 0, circMax);
        ana.setHideOutputImage(true);
        ana.analyze(imgBin);
        for (int i = 0; i < rt.getCounter(); i++) {
           fibArea += rt.getValue("Area", i);
        }
        System.out.println("Fibrosis ="+String.format("%.2f", fibArea)+" µm2");
        // Save masks
        if (saveMask) {
            ImagePlus imgMask = ana.getOutputImage();
            IJ.setRawThreshold(imgMask, 1, 65535, null);
            Prefs.blackBackground = false;
            IJ.run(imgMask, "Convert to Mask", "");
            IJ.run(imgMask, "Create Selection", "");
            Roi roi = imgMask.getRoi();
            Overlay overlay = new Overlay(roi);
            overlay.drawBackgrounds(true);
            img.setOverlay(overlay);
            FileSaver imgSave = new FileSaver(img);
            imgSave.saveAsTiff(path);
            imgMask.close();
        }
        imgBin.close();
        img.close();
        return(fibArea);
    } 
    
  
        
    /**
     * 
     * @param img
     * @param imgSizeX
     * @param path
     * @return 
     */    
    
    public static double findTissueArea(ImagePlus img, double imgScale, boolean zeiss, String path) {
        double pixelSize = img.getCalibration().pixelWidth;
        /* For Zeiss set image calibration to 1 pixel, pyramidal resolution is wrong 
        * it corresponds to high reolution only
        */
        Calibration cal = new Calibration();
        if (zeiss) {
            cal.pixelWidth = 1;
            cal.pixelHeight= 1;
            img.setCalibration(cal);
        }
        
        double tissueArea = 0;
        ResultsTable rt = new ResultsTable();
        IJ.setAutoThreshold(img, "Default");
        Analyzer ana = new Analyzer(img, Measurements.AREA+Measurements.LIMIT, rt);
        ana.measure();
        // rescale area to high resolution image
        if (zeiss)
            tissueArea = rt.getValue("Area", 0) * imgScale * Math.pow(pixelSize, 2);
        else
            tissueArea = rt.getValue("Area", 0);
        System.out.println("Tissue area ="+String.format("%.2f", tissueArea)+" µm2");
        if (zeiss) {
                cal.pixelWidth = cal.pixelHeight = imgScale * Math.pow(pixelSize, 2);
                cal.setUnit("µm");
                img.setCalibration(cal);
        }
        // Save masks
        if (saveMask) {
            IJ.run(img, "Convert to Mask", "");
            FileSaver imgSave = new FileSaver(img);
            imgSave.saveAsTiff(path);
        }
        return(tissueArea);
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

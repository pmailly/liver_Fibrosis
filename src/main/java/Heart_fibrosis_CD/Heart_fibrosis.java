/* Extract red/yellow color define in roi file 
* roi1 = heart
* roi2 = fibrose red
* roi3 = yellow heart no ventricul
* roi4 = background
* Author Philippe Mailly
*/



package Heart_fibrosis_CD;

import static Heart_fibrosis_CD_tools.Heart_fibrosist_tools.analyse;
import static Heart_fibrosis_CD_tools.Heart_fibrosist_tools.findColorMatrix;
import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import static Heart_fibrosis_CD_tools.Heart_fibrosist_tools.find_roi;
import static Heart_fibrosis_CD_tools.Heart_fibrosist_tools.writeHeaders;
import ij.ImageStack;
import ij.gui.Roi;
import static ij.plugin.filter.Analyzer.setOption;
import ij.process.ImageConverter;
import org.apache.commons.io.FilenameUtils;
import Heart_fibrosis_CD_tools.StainMatrix;
import ij.measure.Calibration;
/**
 *
 * @author phm
 */
public class Heart_fibrosis implements PlugIn {
    
    private final boolean canceled = false;
    public String outDirResults, imageDir;
    public BufferedWriter outPutResults;
    

    /**
     * 
     * 
     */
    @Override
    public void run(String arg) {
        try {
            if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            imageDir = IJ.getDirectory("Choose images folder");
            File inDir = new File(imageDir);
            String[] imageFile = inDir.list();
            if (imageFile == null) {
                return;
            }
            
            // create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta;
            meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            int imageNum = 0;
            int series;
            String imageName = "";
            String rootName = "";
            imageDir += File.separator;
            
            
            Arrays.sort(imageFile);
            for (String f : imageFile) {
                String fileExt = FilenameUtils.getExtension(f);
                if (fileExt.equals("czi")) {
                    imageName = imageDir + f;
                    rootName = FilenameUtils.removeExtension(f);
                    // find zip roi file
                    String roiExt =  find_roi(imageFile, rootName);
                    String roiFile = rootName+"#"+ roiExt;
                    if (!new File(imageDir+roiFile).exists()) {
                        System.out.println("No roi file found !!!");
                        return;
                    }
                    
                    series = Integer.parseInt(FilenameUtils.removeExtension(roiExt))-1;
                    imageNum++;
                    reader.setId(imageName);
                    reader.setSeries(series);
                    // get calibration
                    Calibration cal = new Calibration();
                    cal.pixelWidth = meta.getPixelsPhysicalSizeX(series).value().doubleValue();
                    cal.pixelHeight = cal.pixelWidth;
                    cal.pixelDepth = 1;
                    cal.setUnit("microns");
                    System.out.println("x cal = " +cal.pixelWidth+", z cal=" + cal.pixelDepth);

                    if (imageNum == 1) {
                        // create output folder
                        outDirResults = inDir + File.separator+ "Results"+ File.separator;
                        File outDir = new File(outDirResults);
                        if (!Files.exists(Paths.get(outDirResults))) {
                            outDir.mkdir();
                        }
                        /**
                         * Write headers results for results file
                         */

                        String resultsName = "results.xls";
                        String header = "Image\tFibrosis area (red)\tNo fibrosis area (yellow)\tHeart area\n";
                        outPutResults = writeHeaders(outDirResults, resultsName, header);
                    }
                    
                    RoiManager rm = new RoiManager(false);
                    rm.runCommand("Open", imageDir+roiFile);
                    ImporterOptions options = new ImporterOptions();
                    options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                    options.setId(imageName);    
                    options.setQuiet(true);
                    options.setSeriesOn(series, true);
                    // Open images with series
                    ImagePlus img = BF.openImagePlus(options)[0];
                    new ImageConverter(img).convertRGBStackToRGB();
                    img.setCalibration(cal);
                    
                    // Find color matrix with rois
                    StainMatrix stMat = new StainMatrix();
                    double[] matrixValue=new double[9];
                    double [] rgbOD = new double[3];
                    for (int i = 0; i < 3; i++){
                        int r = i +1;
                        findColorMatrix(img, rgbOD, rm.getRoi(r));
                        matrixValue[3*i+0]=rgbOD[0];
                        matrixValue[3*i+1]=rgbOD[1];
                        matrixValue[3*i+2]=rgbOD[2];
                    }
                    stMat.init("From ROI", matrixValue[0], matrixValue[1], matrixValue[2],
                              matrixValue[3], matrixValue[4], matrixValue[5],
                              matrixValue[6], matrixValue[7], matrixValue[8]);
                    
                    // Do color deconvolution
                    ImageStack[] imgStack = stMat.compute(false, true, img);
                    img.close();
                    
                    // select Red channel
                    ImagePlus imgRed = new ImagePlus("Fibrosis",imgStack[0]);
                    imgRed.setCalibration(cal);
                    IJ.run(imgRed, "Median...", "radius=4 stack");
                    IJ.setAutoThreshold(imgRed, "Otsu");
                    setOption("BlackBackground", false);
                    Roi roi = rm.getRoi(0);
                    imgRed.setRoi(roi);
                    imgRed.updateAndDraw();
                    double heartArea = roi.getStatistics().area;
                    System.out.println("heart ="+heartArea);
                    double redArea = analyse(imgRed, outDirResults+rootName+"_Fibrosis_Mask.tif");

                    // select Yellow channel
                    ImagePlus imgYellow = new ImagePlus("Heart",imgStack[1]);
                    imgYellow.setCalibration(cal);
                    IJ.run(imgYellow, "Median...", "radius=4 stack");
                    IJ.setAutoThreshold(imgYellow, "Otsu");
                    setOption("BlackBackground", false);
                    imgYellow.setRoi(roi); 
                    imgYellow.updateAndDraw();
                    
                    double yellowArea = analyse(imgYellow, outDirResults+rootName+"_Heart_Mask.tif");
                    outPutResults.write(rootName+"\t"+redArea+"\t"+ yellowArea + "\t" + heartArea+ "\n");
                    outPutResults.flush();
                }
            }
            if (Files.exists(Paths.get(outDirResults+"results.xls")))
                outPutResults.close();
                

        } catch (IOException | DependencyException | ServiceException | FormatException ex) {
            Logger.getLogger(Heart_fibrosis.class.getName()).log(Level.SEVERE, null, ex);
        }

        IJ.showStatus("Process done ...");

    }
    
}

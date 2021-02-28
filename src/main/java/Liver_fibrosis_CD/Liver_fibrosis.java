/* Extract brown color define in matrix color 
* from tunnel image coloration
* Author Philippe Mailly
*/



package Liver_fibrosis_CD;

import static Liver_fibrosis_CD_tools.Liver_fibrosis_tools.clij_analyse;
import static Liver_fibrosis_CD_tools.Liver_fibrosis_tools.color_meth;
import static Liver_fibrosis_CD_tools.Liver_fibrosis_tools.dialog;
import static Liver_fibrosis_CD_tools.Liver_fibrosis_tools.findTissueArea;
import static Liver_fibrosis_CD_tools.Liver_fibrosis_tools.tileSize;
import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;
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
import static Liver_fibrosis_CD_tools.Liver_fibrosis_tools.writeHeaders;
import Liver_fibrosis_CD_tools.StainMatrix;
import ij.ImageStack;
import ij.process.ImageConverter;
import org.apache.commons.io.FilenameUtils;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.process.ImageStatistics;
import loci.common.Region;





/**
 *
 * @author phm
 */
public class Liver_fibrosis implements PlugIn {
    
    private final boolean canceled = false;
    public String outDirResults, imageDir;
    public BufferedWriter outPutResults;
    // Red Sirius Input Vector Matrix
    private final double[] MODx = new double[3];
    private final double[] MODy = new double[3];
    private final double[] MODz = new double[3];
    private boolean zeiss = false;

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
            imageDir = dialog();
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
            int series = 0;
            String imageName = "";
            String rootName = "";
            imageDir += File.separator;
            Arrays.sort(imageFile);
            for (String f : imageFile) {
                String fileExt = FilenameUtils.getExtension(f);
                if (fileExt.equals("ndpi") || fileExt.equals("czi")) {
                    if (fileExt.equals("czi"))
                        zeiss =true;
                    imageName = imageDir + f;
                    rootName = FilenameUtils.removeExtension(f);
                    imageNum++;
                    reader.setId(imageName);
                    reader.setSeries(series);
                    // get calibration
                    Calibration cal = new Calibration();
                    cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
                    cal.pixelHeight = cal.pixelWidth;
                    cal.pixelDepth = 1;
                    cal.setUnit(meta.getPixelsPhysicalSizeX(0).unit().getSymbol());
                    System.out.println("x cal = " +cal.pixelWidth+", "+cal.getUnit());

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
                        String header = "Image\tFibrosis area\tTissue area\t%Fibrosis\n";
                        outPutResults = writeHeaders(outDirResults, resultsName, header);
                        
                        // Initialze color matrix for H-DAB
                        MODx[0] = 0.65;
			MODy[0] = 0.704;
			MODz[0] = 0.286;
 
			MODx[1] = 0.268;
			MODy[1] = 0.57;
			MODz[1] = 0.776;
 
			MODx[2] = 0.0;
			MODy[2] = 0.0;
			MODz[2] = 0.0;


                    }
                    int imgSizeX = reader.getSizeX();
                    int imgSizeY = reader.getSizeY();
//                    int tileSizeX = reader.getOptimalTileWidth();
//                    int tileSizeY = reader.getOptimalTileHeight();
                    int tileSizeX = (int)tileSize;
                    int tileSizeY = (int)tileSize;
                    int nXTiles = imgSizeX / tileSizeX;
                    int nYTiles = imgSizeY / tileSizeY;
                    System.out.println("Tile size X = "+tileSizeX+" tile size Y = "+tileSizeY);
                    
                    // Find tissue area on last series (lowest resolution)
                    int nbSeries = reader.getSeriesCount();
                    
                    /* Case Zeiss images the calibration if only correct
                    * for series 0. Resolution facteur for other series  2^series
                    */
                    
                    //get lowest pyramidal resolution
                    series = nbSeries - 3;
                    reader.setSeries(series);
                    ImporterOptions options = new ImporterOptions();
                    options.setColorMode(ImporterOptions.COLOR_MODE_DEFAULT);
                    options.setId(imageName);    
                    options.setQuiet(true);
                    options.setCrop(true);
                    options.setSeriesOn(series, true);
                    
                    // remove first line for czi (black)
                    int xline = 0;
                    if (zeiss)
                        xline = 1;
                    Region reg = new Region(0, xline, reader.getSizeX() - 1, reader.getSizeY() - 2);
                    options.setCropRegion(series, reg);
                    options.doCrop();
                    
                    // Open images
                    ImagePlus imgLow = BF.openImagePlus(options)[0];
                    new ImageConverter(imgLow).convertToRGB();
                    new ImageConverter(imgLow).convertToGray16();
                    // For Zeiss image image pyramidal calibration is wrong
                    if (zeiss) {
                        double lowPixelSizeX = cal.pixelWidth * Math.pow(2, nbSeries - 3);
                        double lowPixelSizeY = lowPixelSizeX;
                        Calibration lowCal = new Calibration();
                        lowCal.pixelWidth = lowPixelSizeX;
                        lowCal.pixelHeight = lowPixelSizeY;
                        imgLow.setCalibration(lowCal);
                    }
                    double tissueArea = findTissueArea(imgLow, outDirResults+rootName+"_Tissue_Mask.tif");
                    imgLow.close();

                    // Find tiles dimensions
                    int currentTile = 0;
                    double fibrosisArea = 0;
                    series = 0;
                    reader.setSeries(series);
                    for (int y = 0; y < nYTiles; y++) {
                        for (int x = 0; x <nXTiles; x++) {
                            currentTile++;
                            IJ.showStatus("opening tile "+ currentTile + "/"+nYTiles+nXTiles);
                        // The x and y coordinates for the current tile
                            int tileX = x * tileSizeX;
                            int tileY = y * tileSizeY;
                            reg = new Region(tileX, tileY, tileSizeX, tileSizeY);
                            options.setSeriesOn(series, true);
                            options.setCropRegion(series, reg);
                            options.doCrop();
                            // Open images
                            ImagePlus img = BF.openImagePlus(options)[0];
                            
                            /* if no object do nothing 
                            * for Zeiss black image if no object
                            * for Hamamatsu scan all part
                            */
                            ImageStatistics stats = img.getStatistics(ImageStatistics.MEAN+ImageStatistics.STD_DEV);
                            if ((zeiss && stats.mean >= 10) || (!zeiss && stats.stdDev > 10)) {
                                new ImageConverter(img).convertToRGB();
                                ImagePlus imgFibrosis = new ImagePlus();
                                if (color_meth.equals("Tunel")) {
                                  StainMatrix stMat = new StainMatrix();
                                  stMat.init("From Value", MODx[0], MODy[0], MODz[0],
                                  MODx[1], MODy[1], MODz[1], MODx[2], MODy[2], MODz[2]);
                                  ImageStack[] imgStack = stMat.compute(false, true, img);
                                  imgFibrosis = new Duplicator().run(new ImagePlus("Fibrosis",imgStack[1]));
                                }
                                else {
                                    imgFibrosis = new Duplicator().run(img);
                                    new ImageConverter(imgFibrosis).convertToGray16();
                                }
                                imgFibrosis.setCalibration(cal);
                                
                                // Find fibrosis area
                                if (color_meth.equals("Tunel")) {
                                    if (zeiss)
                                        fibrosisArea += clij_analyse(img, imgFibrosis, "IsoData", "log", 0, 1000, 500000, outDirResults+rootName+"_Fibrosis_Mask_"+tileX+"_"+ tileY+".tif");
                                    else
                                        fibrosisArea += clij_analyse(img, imgFibrosis, "MaxEntropy", "exp", 6, 2000, Double.MAX_VALUE, outDirResults+rootName+"_Fibrosis_Mask_"+tileX+"_"+ tileY+".tif");
                                }
                                else
                                    fibrosisArea +=clij_analyse(img, imgFibrosis, "Triangle", "exp", 4, 10, 100000, outDirResults+rootName+"_Fibrosis_Mask_"+tileX+"_"+ tileY+".tif");
                                imgFibrosis.close();
                            }
                            else
                                img.close();
                        }
                    }
                    outPutResults.write(rootName+"\t"+fibrosisArea+"\t"+ tissueArea + "\t"+ (fibrosisArea / tissueArea)*100+"\n");
                    outPutResults.flush();
                }
            }
            if (Files.exists(Paths.get(outDirResults+"results.xls")))
                outPutResults.close();
                

        } catch (IOException | DependencyException | ServiceException | FormatException ex) {
            Logger.getLogger(Liver_fibrosis.class.getName()).log(Level.SEVERE, null, ex);
        }

        IJ.showStatus("Process done ...");

    }
    
}

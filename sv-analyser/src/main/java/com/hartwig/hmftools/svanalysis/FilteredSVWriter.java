package com.hartwig.hmftools.svanalysis;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class FilteredSVWriter {

    private final String mOutputPath;
    private final String mVcfFileLocation;
    private BufferedWriter mFileWriter;

    private boolean mRunPONFilter;
    private boolean mLogInsSVs;

    private static final Logger LOGGER = LogManager.getLogger(FilteredSVWriter.class);

    public FilteredSVWriter(final String vcfFileLocation, final String outputPath)
    {
        mRunPONFilter = false;
        mLogInsSVs = false;

        mOutputPath = outputPath;
        mVcfFileLocation = vcfFileLocation;
        mFileWriter = null;
    }

    public void setRunPONFilter(boolean toggle) { mRunPONFilter = toggle; }
    public void setLogInsSVs(boolean toggle) { mLogInsSVs = toggle; }

    public void processVcfFiles()
    {
        final List<File> vcfFiles;

        final Path root = Paths.get(mVcfFileLocation);

        try (final Stream<Path> stream = Files.walk(root, 5, FileVisitOption.FOLLOW_LINKS)) {

            vcfFiles = stream.map(p -> p.toFile())
                    .filter(p -> !p.isDirectory())
                    .filter(p_-> p_.getName().endsWith("somaticSV_bpi.vcf"))
                    .collect(Collectors.toList());


            LOGGER.info("found {} BPI VCF files", vcfFiles.size());

                // add the filtered and passed SV entries for each file
            for(final File vcfFile : vcfFiles)
            {
                if(vcfFile.isDirectory())
                    continue;

                if(!vcfFile.getPath().contains("structuralVariants/bpi/"))
                    continue;

                if(!vcfFile.getName().endsWith("somaticSV_bpi.vcf"))
                    continue;

                LOGGER.debug("BPI VCF path({}) file({})", vcfFile.getPath(), vcfFile.getName());

                // extract sampleId from the directory or file name
                String[] itemsStr = vcfFile.getName().split("_");

                if(itemsStr.length != 4)
                continue;

                String sampleId = itemsStr[1];
                LOGGER.debug("sampleId({})", sampleId);

                List<StructuralVariant> variants = readFromVcf(vcfFile.getPath());

                if(mRunPONFilter)
                    generateFilteredSVFile(variants, sampleId);

                if(mLogInsSVs)
                    logInsData(variants, sampleId);
            }

            if(mFileWriter != null)
                mFileWriter.close();

        }
        catch(Exception e)
        {

        }
    }

    private List<StructuralVariant> readFromVcf(String vcfFilename)  {
        final StructuralVariantFactory factory = new StructuralVariantFactory(false);
        try (final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(vcfFilename,
                new VCFCodec(),
                false)) {
            reader.iterator().forEach(factory::addVariantContext);
        }
        catch(Exception e) {
        }

        return factory.results();
    }


    private void generateFilteredSVFile(final List<StructuralVariant> variants, final String sampleId)
    {
        try {
            if(mFileWriter == null)
            {
                String outputFileName = mOutputPath;
                if(!outputFileName.endsWith("/"))
                    outputFileName += "/";

                outputFileName += "svs_incl_filtered.csv";

                Path outputFile = Paths.get(outputFileName);

                mFileWriter = Files.newBufferedWriter(outputFile);

                mFileWriter.write("SampleId,SvId,Type,ChrStart,PosStart,OrientStart,ChrEnd,PosEnd,OrientEnd,Filters\n");
            }

            for(final StructuralVariant var : variants)
            {
                String filtersStr = var.filter();

                if(filtersStr.equals("PASS") || filtersStr.equals("[]") || filtersStr.equals(".") || filtersStr.isEmpty())
                {
                    LOGGER.debug("var({}) was a PASS", var.id());
                    filtersStr = "PASS";
                }
                else
                {
                    // make tokenisable for further searches
                    LOGGER.debug("var({}) has filters: {}", var.id(), var.filter());

                    if(filtersStr.charAt(0) == '[')
                        filtersStr = filtersStr.substring(1);
                    if(filtersStr.charAt(filtersStr.length() - 1) == ']')
                        filtersStr = filtersStr.substring(0, filtersStr.length() - 1);
                    if(!filtersStr.isEmpty())
                        filtersStr = filtersStr.replace(",", ";");
                }

                mFileWriter.write(
                        String.format("%s,%s,%s,%s,%d,%d,%s,%d,%d,%s",
                                sampleId, var.id(), var.type(),
                                var.chromosome(true), var.position(true), var.orientation(true),
                                var.chromosome(false), var.position(false), var.orientation(false), filtersStr));

                mFileWriter.newLine();
            }

        }
        catch (final IOException e) {
            LOGGER.error("error writing to outputFile");
        }
    }

    private void logInsData(List<StructuralVariant> variants, final String sampleId)
    {
        try {

            if(mFileWriter == null)
            {
                String outputFileName = mOutputPath;
                if(!outputFileName.endsWith("/"))
                    outputFileName += "/";

                outputFileName += "sv_inserts.csv";

                Path outputFile = Paths.get(outputFileName);

                mFileWriter = Files.newBufferedWriter(outputFile);

                mFileWriter.write("SampleId,SvId,Type,ChrStart,PosStart,OrientStart,ChrEnd,PosEnd,OrientEnd,InsertSeq\n");
            }

            for(final StructuralVariant var : variants)
            {
                if(var.type() != StructuralVariantType.INS)
                    continue;

                String filtersStr = var.filter();

                if(!filtersStr.equals("PASS") && !filtersStr.equals(".") && !filtersStr.equals("[]"))
                {
                    continue;
                }

                mFileWriter.write(
                        String.format("%s,%s,%s,%s,%d,%d,%s,%d,%d,%s",
                                sampleId, var.id(), var.type(),
                                var.chromosome(true), var.position(true), var.orientation(true),
                                var.chromosome(false), var.position(false), var.orientation(false), var.insertSequence()));

                mFileWriter.newLine();
            }

        }
        catch (final IOException e) {
            LOGGER.error("error writing to outputFile");
        }
    }
}

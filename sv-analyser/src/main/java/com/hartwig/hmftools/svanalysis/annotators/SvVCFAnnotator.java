package com.hartwig.hmftools.svanalysis.annotators;

import java.io.BufferedWriter;
import java.io.File;
import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class SvVCFAnnotator {

    private final String mOutputPath;
    private final String mVcfFileLocation;
    private BufferedWriter mFileWriter;

    private static final Logger LOGGER = LogManager.getLogger(SvVCFAnnotator.class);

    public SvVCFAnnotator(final String vcfFileLocation, final String outputPath)
    {
        mOutputPath = outputPath;
        mVcfFileLocation = vcfFileLocation;
        mFileWriter = null;
    }

    public void processVcfFiles()
    {
        final List<File> vcfFiles;

        final Path root = Paths.get(mVcfFileLocation);

        try (final Stream<Path> stream = Files.walk(root, 5, FileVisitOption.FOLLOW_LINKS)) {

            vcfFiles = stream.map(p -> p.toFile())
                    .filter(p -> !p.isDirectory())
                    .filter(p_-> p_.getName().endsWith("somaticSV_bpi.vcf"))
                    .collect(Collectors.toList());


            LOGGER.debug("found {} BPI VCF files", vcfFiles.size());

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

                List<VariantContext> variants = readFromVcf(vcfFile.getPath());

                for(final VariantContext variant : variants)
                {
                    if(variant.hasAttribute("IMPRECISE"))
                    {
                        LOGGER.debug("var({}) has imprecise call", variant.getID());
                    }
                }
                // generateFilteredSVFile(variants, sampleId);
            }

            if(mFileWriter != null)
                mFileWriter.close();

        }
        catch(Exception e)
        {

        }
    }

    private List<VariantContext> readFromVcf(String vcfFilename)
    {
        List<VariantContext> variants = Lists.newArrayList();

        final StructuralVariantFactory factory = new StructuralVariantFactory(false);
        try (final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(vcfFilename,
                new VCFCodec(),
                false))
        {
            while(reader.iterator().hasNext())
                variants.add(reader.iterator().next());
        }
        catch(Exception e) {
        }

        return variants;
    }

}

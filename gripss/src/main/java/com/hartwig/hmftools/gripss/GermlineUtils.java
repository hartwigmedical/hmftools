package com.hartwig.hmftools.gripss;

import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class GermlineUtils
{
    public static final Logger GM_LOGGER = LogManager.getLogger(GripssApplication.class);

    public static String stripBam(final String sampleId)
    {
        return sampleId.replaceAll("_dedup.realigned.bam","")
                .replaceAll(".sorted", "")
                .replaceAll(".bam", "");
    }

    public static List<String> findVcfFiles(final String batchRunRootDir)
    {
        // current prod examples
        // structuralVariants/gridss/CPCT02030278R_CPCT02030278T/CPCT02030278R_CPCT02030278T.gridss.vcf.gz
        // structural_caller/WIDE01010356T.gridss.unfiltered.vcf.gz
        final List<String> vcfFiles = Lists.newArrayList();

        try
        {
            final Stream<Path> stream = Files.walk(Paths.get(batchRunRootDir), 5, FileVisitOption.FOLLOW_LINKS);

            vcfFiles.addAll(stream.filter(x -> !x.toFile().isDirectory())
                    .map(x -> x.toFile().toString())
                    .filter(x -> matchesGridssVcf(x))
                    .collect(Collectors.toList()));

            GM_LOGGER.info("found {} VCF files", vcfFiles.size());
        }
        catch (Exception e)
        {
            GM_LOGGER.error("failed find directories for batchDir({}) run: {}", batchRunRootDir, e.toString());
        }

        return vcfFiles;
    }

    private static boolean matchesGridssVcf(final String filename)
    {
        return filename.endsWith(".gridss.vcf") || filename.endsWith(".gridss.unfiltered.vcf")
                || filename.endsWith(".gridss.vcf.gz") || filename.endsWith(".gridss.unfiltered.vcf.gz");
    }

}

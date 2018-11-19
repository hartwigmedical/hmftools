package com.hartwig.hmftools.bachelorpp;

import java.io.File;
import java.io.IOException;
import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bachelorpp.types.BachelorGermlineVariant;
import com.hartwig.hmftools.common.pileup.Pileup;
import com.hartwig.hmftools.common.pileup.PileupFile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class AlleleDepthLoader {

    private static final String MPILEUP_FILE_EXTN = ".mpu";

    private static final Logger LOGGER = LogManager.getLogger(AlleleDepthLoader.class);

    private String mSampleId;
    private List<Pileup> mPileupData;

    public AlleleDepthLoader()
    {
        mSampleId = "";
        mPileupData = Lists.newArrayList();
    }

    public void setSampleId(final String sampleId) {
        mSampleId = sampleId;
    }

    public boolean loadMiniPileupData(final String mpDirectory)
    {
        final List<File> mpuFiles;

        final Path root = Paths.get(mpDirectory);

        try (final Stream<Path> stream = Files.walk(root, 1, FileVisitOption.FOLLOW_LINKS))
        {
            mpuFiles = stream.map(p -> p.toFile())
                    .filter(p -> !p.isDirectory())
                    .filter(p_ -> p_.getName().endsWith(MPILEUP_FILE_EXTN))
                    .collect(Collectors.toList());

            LOGGER.debug("found {} mini-pileup files", mpuFiles.size());

            for (final File mpuFile : mpuFiles) {

                if (mpuFile.isDirectory()) {
                    continue;
                }

                LOGGER.debug("processing mini-pileup file({})", mpuFile.getName());

                List<Pileup> pileups = PileupFile.read(mpuFile.getPath());
                mPileupData.addAll(pileups);
            }

            LOGGER.info("loaded {} pileup files, {} records", mpuFiles.size(), mPileupData.size());
        }
        catch (IOException e)
        {
            LOGGER.error("failed to process pileup files from dir({})", mpDirectory);
            return false;
        }

        return true;
    }

    public boolean applyPileupData(List<BachelorGermlineVariant> germlineVariants)
    {
        // match up read info from MP with the bachelor records
        for (BachelorGermlineVariant bachRecord : germlineVariants)
        {
            if(bachRecord.isReadDataSet())
                continue;

            boolean matched = false;

            for (final Pileup pileup : mPileupData)
            {
                if (!bachRecord.chromosome().equals(pileup.chromosome()))
                    continue;

                if (bachRecord.position() != pileup.position())
                    continue;

                matched = true;
                bachRecord.setRefCount(pileup.referenceCount());

                if (pileup.insertCount() > 0)
                {
                    bachRecord.setAltCount(pileup.insertCount());
                }
                else if (pileup.deleteCount() > 0)
                {
                    bachRecord.setAltCount(pileup.deleteCount());
                }
                else if (bachRecord.alts().length() == bachRecord.ref().length() && bachRecord.alts().length() == 1)
                {
                    if (bachRecord.alts().charAt(0) == 'A')
                    {
                        bachRecord.setAltCount(pileup.aMismatchCount());
                    }
                    else if (bachRecord.alts().charAt(0) == 'C')
                    {
                        bachRecord.setAltCount(pileup.cMismatchCount());
                    }
                    else if (bachRecord.alts().charAt(0) == 'G')
                    {
                        bachRecord.setAltCount(pileup.gMismatchCount());
                    }
                    else if (bachRecord.alts().charAt(0) == 'T')
                    {
                        bachRecord.setAltCount(pileup.tMismatchCount());
                    }
                }

                LOGGER.debug("sample({} chr({}) position({}) matched, counts(ref={} alt={})",
                        mSampleId, bachRecord.chromosome(), bachRecord.position(), bachRecord.getRefCount(), bachRecord.getAltCount());
            }

            if(!matched)
            {
                LOGGER.warn("sample({} var({}) chr({}) position({}) no pile-up record found",
                        mSampleId, bachRecord.variantId(), bachRecord.chromosome(), bachRecord.position());
                return false;
            }
        }

        return true;
    }

    public final List<Pileup> getPileupData() {
        return mPileupData;
    }

}

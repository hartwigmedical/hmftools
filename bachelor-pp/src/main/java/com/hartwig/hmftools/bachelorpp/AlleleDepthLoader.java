package com.hartwig.hmftools.bachelorpp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
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

    private static final int BACHELOR_CSV_FIELD_COUNT = 15;
    private static final String MPILEUP_FILE_EXTN = ".mpu";

    private static final Logger LOGGER = LogManager.getLogger(AlleleDepthLoader.class);

    private String sampleId;
    private List<BachelorGermlineVariant> bachelorVariants;
    private List<Pileup> pileupData;

    public AlleleDepthLoader()
    {
        sampleId = "";
        bachelorVariants = Lists.newArrayList();
        pileupData = Lists.newArrayList();
    }

    public void setSampleId(final String sampleId) {
        this.sampleId = sampleId;
    }

    public boolean loadBachelorMatchData(final String filename)
    {
        if (filename.isEmpty())
            return false;

        try {

            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine(); // skip header

            while ((line = fileReader.readLine()) != null) {

                // parse CSV data
                String[] items = line.split(",");

                // CSV fields PATIENT,SOURCE,PROGRAM,ID,GENES,TRANSCRIPT_ID,CHROM,POS,REF,ALTS,EFFECTS

                if (items.length != BACHELOR_CSV_FIELD_COUNT) {
                    LOGGER.warn("invalid item count({}), recordIndex({}) in file({})", items.length, bachelorVariants.size(), filename);
                    return false;
                }

                final String patientId = items[0];

                if (!sampleId.equals("*") && !sampleId.equals("") && !sampleId.contains(patientId)) {
                    continue;
                }

                BachelorGermlineVariant bachRecord = new BachelorGermlineVariant(patientId,
                        items[1],
                        items[2],
                        items[3],
                        items[4],
                        items[5],
                        items[6],
                        Long.parseLong(items[7]),
                        items[8],
                        items[9],
                        items[10],
                        items[11],
                        items[12],
                        Boolean.parseBoolean(items[13]),
                        Integer.parseInt(items[14]));

                bachelorVariants.add(bachRecord);
            }

            LOGGER.debug("loaded {} bachelor records", bachelorVariants.size());

        }
        catch (IOException exception)
        {
            LOGGER.error("Failed to read bachelor input CSV file({})", filename);
            return false;
        }

        return true;
    }

    public final List<BachelorGermlineVariant> getBachelorVariants() {
        return bachelorVariants;
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
                pileupData.addAll(pileups);
            }

            LOGGER.info("loaded {} pileup files, {} records", mpuFiles.size(), pileupData.size());
        }
        catch (IOException e)
        {
            LOGGER.error("failed to process pileup files from dir({})", mpDirectory);
            return false;
        }

        applyPileupData();

        return true;
    }

    private void applyPileupData()
    {
        // match up read info from MP with the bachelor records
        for (BachelorGermlineVariant bachRecord : bachelorVariants)
        {
            for (final Pileup pileup : pileupData)
            {
                if (!bachRecord.chromosome().equals(pileup.chromosome()))
                    continue;

                if (bachRecord.position() != pileup.position())
                    continue;

                bachRecord.setRefCount(pileup.referenceCount());

                if (pileup.insertions() > 0)
                {
                    bachRecord.setAltCount(pileup.insertions());
                }
                else if (pileup.deletions() > 0)
                {
                    bachRecord.setAltCount(pileup.deletions());
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
                        sampleId, bachRecord.chromosome(), bachRecord.position(), bachRecord.getRefCount(), bachRecord.getAltCount());
            }
        }
    }

    public final List<Pileup> getPileupData() {
        return pileupData;
    }

}

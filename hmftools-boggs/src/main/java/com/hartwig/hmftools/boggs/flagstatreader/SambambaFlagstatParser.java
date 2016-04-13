package com.hartwig.hmftools.boggs.flagstatreader;

import org.jetbrains.annotations.NotNull;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class SambambaFlagstatParser implements FlagstatParser {

    public SambambaFlagstatParser() {
    }

    @NotNull
    public FlagstatData parse(@NotNull File file) throws IOException {
        BufferedReader fileReader  = new BufferedReader(new FileReader(file));
        FlagStatsBuilder qcPassedReadsBuilder = new FlagStatsBuilder();
        FlagStatsBuilder qcFailedReadsBuilder = new FlagStatsBuilder();

        int lineNr = 0;
        String line;
        while ((line = fileReader.readLine()) != null) {
            lineNr++;
            // KODU: Assume format <qcPassedNr> + <qcFailedNr> <description sentence>
            String[] parts = line.split(" ");
            long qcPassed = Long.parseLong(parts[0]);
            long qcFailed = Long.parseLong(parts[2]);

            switch (lineNr) {
                case 1: {
                    qcPassedReadsBuilder.setTotal(qcPassed);
                    qcFailedReadsBuilder.setTotal(qcFailed);
                    break;
                }
                case 2: {
                    qcPassedReadsBuilder.setSecondary(qcPassed);
                    qcFailedReadsBuilder.setSecondary(qcFailed);
                    break;
                }
                case 3: {
                    qcPassedReadsBuilder.setSupplementary(qcPassed);
                    qcFailedReadsBuilder.setSupplementary(qcFailed);
                    break;
                }
                case 4: {
                    qcPassedReadsBuilder.setDuplicates(qcPassed);
                    qcFailedReadsBuilder.setDuplicates(qcFailed);
                    break;
                }
                case 5: {
                    qcPassedReadsBuilder.setMapped(qcPassed);
                    qcFailedReadsBuilder.setMapped(qcFailed);
                    break;
                }
                case 6: {
                    qcPassedReadsBuilder.setPairedInSequencing(qcPassed);
                    qcFailedReadsBuilder.setPairedInSequencing(qcFailed);
                    break;
                }
                case 7: {
                    qcPassedReadsBuilder.setRead1(qcPassed);
                    qcFailedReadsBuilder.setRead1(qcFailed);
                    break;
                }
                case 8: {
                    qcPassedReadsBuilder.setRead2(qcPassed);
                    qcFailedReadsBuilder.setRead2(qcFailed);
                    break;
                }
                case 9: {
                    qcPassedReadsBuilder.setProperlyPaired(qcPassed);
                    qcFailedReadsBuilder.setProperlyPaired(qcFailed);
                    break;
                }
                case 10: {
                    qcPassedReadsBuilder.setItselfAndMateMapped(qcPassed);
                    qcFailedReadsBuilder.setItselfAndMateMapped(qcFailed);
                    break;
                }
                case 11: {
                    qcPassedReadsBuilder.setSingletons(qcPassed);
                    qcFailedReadsBuilder.setSingletons(qcFailed);
                    break;
                }
                case 12: {
                    qcPassedReadsBuilder.setMateMappedToDifferentChr(qcPassed);
                    qcFailedReadsBuilder.setMateMappedToDifferentChr(qcFailed);
                    break;
                }
                case 13: {
                    qcPassedReadsBuilder.setMateMappedToDifferentChrMapQ5(qcPassed);
                    qcFailedReadsBuilder.setMateMappedToDifferentChrMapQ5(qcFailed);
                    break;
                }
            }
        }

        fileReader.close();

        return new FlagstatData(file.getPath(), qcPassedReadsBuilder.build(), qcFailedReadsBuilder.build());
    }
}

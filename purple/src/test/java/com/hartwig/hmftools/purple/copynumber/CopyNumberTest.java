package com.hartwig.hmftools.purple.copynumber;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_2;
import static com.hartwig.hmftools.purple.copynumber.PurpleCopyNumberFactory.calculateDeletedDepthWindows;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleTestUtils;
import com.hartwig.hmftools.common.test.GeneTestUtils;
import com.hartwig.hmftools.purple.gene.GeneCopyNumberBuilder;

import org.junit.Test;

public class CopyNumberTest
{
    @Test
    public void testDeletedDepthWindowsCalc()
    {
        List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();
        copyNumbers.add(createCopyNumber("1", 1, 10000, 2, 10));
        copyNumbers.add(createCopyNumber("1", 10001, 20000, 0.4, 10));
        copyNumbers.add(createCopyNumber("1", 20001, 30000, 2, 10));
        copyNumbers.add(createCopyNumber("1", 30001, 40000, 0.1, 10));
        copyNumbers.add(createCopyNumber("1", 40001, 50000, 2, 60));

        copyNumbers.add(createCopyNumber("9", 8000000, 10000000, 0.4, 40));
        copyNumbers.add(createCopyNumber("9", 10000001, 11000000, 2, 40));
        copyNumbers.add(createCopyNumber("9", 11000001, 13000000, 0.4, 20));

        // total = 200

        double deletedPercent = calculateDeletedDepthWindows(copyNumbers);
        assertEquals(0.25, deletedPercent, 0.01);

    }
    private static PurpleCopyNumber createCopyNumber(String chromosome, int start, int end, double copyNumber, int depthWindowCount)
    {
        return PurpleTestUtils.createCopyNumber(chromosome, start, end, copyNumber)
                .depthWindowCount(depthWindowCount)
                .build();
    }
}

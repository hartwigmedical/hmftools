package com.hartwig.hmftools.pavereverse.batch;

import static com.hartwig.hmftools.pavereverse.batch.BaseSequenceVariantsEncoder.CHROMOSOME;
import static com.hartwig.hmftools.pavereverse.batch.BaseSequenceVariantsEncoder.TRANSCRIPT;
import static com.hartwig.hmftools.pavereverse.batch.BaseSequenceVariantsEncoder.VARIANTS;
import static com.hartwig.hmftools.pavereverse.batch.BatchSequenceVariantsEncoderTest.parseChanges;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Objects;
import java.util.Set;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.pavereverse.BaseSequenceChange;
import com.hartwig.hmftools.pavereverse.ReversePaveTestBase;

import org.apache.commons.io.FileUtils;
import org.junit.Assert;
import org.junit.Test;

public class BatchProcessorTest extends ReversePaveTestBase
{
    @Test
    public void processTest()
    {
        ClassLoader classLoader = getClass().getClassLoader();
        File input = new File(Objects.requireNonNull(classLoader.getResource("batch/test1.tsv")).getFile());
        File outputDir = new File("target/test-output/batch/");
        File output = new File(outputDir, "test1.tsv");
        new BatchProcessor(reversePave).process(input.getAbsolutePath(), output.getAbsolutePath());

        try(DelimFileReader outputReader = new DelimFileReader(output.getAbsolutePath()))
        {
            List<DelimFileReader.Row> resultRows = outputReader.stream().toList();
            assertEquals(4, resultRows.size());

            // First is MTOR L2230V, see ReversePaveTest
            assertEquals("ENST00000361445", resultRows.get(0).get(TRANSCRIPT));
            assertEquals("1", resultRows.get(0).get(CHROMOSOME));
            String variantsStr = resultRows.get(0).get(VARIANTS);
            Set<BaseSequenceChange> parsedChanges = parseChanges(variantsStr, "1");
            assertEquals(4, parsedChanges.size());
            Assert.assertTrue(parsedChanges.contains(new BaseSequenceChange("A", "C", "1", 11_122_101)));

            // 2nd had no transcript, was BRAF V600E
            assertEquals("ENST00000646891", resultRows.get(1).get(TRANSCRIPT));
            assertEquals("7", resultRows.get(1).get(CHROMOSOME));
            variantsStr = resultRows.get(1).get(VARIANTS);
            parsedChanges = parseChanges(variantsStr, "7");
            assertEquals(2, parsedChanges.size());
            Assert.assertTrue(parsedChanges.contains(new BaseSequenceChange("A", "T", "7", 140753336)));
        }
        FileUtils.deleteQuietly(output);
    }
}

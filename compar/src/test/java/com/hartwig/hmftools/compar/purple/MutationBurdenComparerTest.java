package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.compar.common.CategoryType.MUTATION_BURDEN;
import static com.hartwig.hmftools.compar.common.MismatchType.VALUE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.common.ExpectedCache;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.SourceData;
import com.hartwig.hmftools.compar.common.SourceType;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

public class MutationBurdenComparerTest
{
    private static final String SAMPLE_ID_1 = "CLAUDETESTSAMPLE1";
    private static final String SAMPLE_ID_2 = "CLAUDETESTSAMPLE2";

    @Rule
    public TemporaryFolder TempDir = new TemporaryFolder();

    @Test
    public void mutationBurdenComparedApartFromPurity() throws IOException
    {
        File expectedFile = TempDir.newFile();
        Files.write(
                expectedFile.toPath(), List.of(
                        String.join(TSV_DELIM, "SampleId", "Category", "Key", "FieldValues"),
                        String.join(TSV_DELIM, SAMPLE_ID_1, "PURITY", "PURITY", "Purity=0.8"),
                        String.join(TSV_DELIM, SAMPLE_ID_1, "MUTATION_BURDEN", "MUTATION_BURDEN", "TmbPerMb=20;MsStatus=MSS"),
                        String.join(TSV_DELIM, SAMPLE_ID_2, "PURITY", "PURITY", "Purity=0.8;TmbPerMb=20")));

        ExpectedCache expectedCache = new ExpectedCache();
        expectedCache.loadFiles(expectedFile.getAbsolutePath());

        File purpleDir = TempDir.newFolder();
        PurityData purity = TestPurityDataBuilder.BUILDER.create();
        PurityContextFile.write(purpleDir.getAbsolutePath(), SAMPLE_ID_1, purity.Purity);
        PurityContextFile.write(purpleDir.getAbsolutePath(), SAMPLE_ID_2, purity.Purity);

        ConfigBuilder configBuilder = new ConfigBuilder();
        PipelineSourcePaths.registerConfig(configBuilder);
        configBuilder.setValue(PURPLE_DIR_CFG + "_" + SourceType.NEW.configStr(), purpleDir.getAbsolutePath());
        PipelineSourcePaths purplePaths = PipelineSourcePaths.fromConfig(SourceType.NEW, configBuilder);

        ComparConfig config = new ComparConfig();
        config.Sources.set(0, SourceData.fromExpectedSource(SourceType.OLD, expectedCache, null));
        config.Sources.set(1, SourceData.fromPipelineSource(SourceType.NEW, purplePaths, null));

        PurityComparer purityComparer = new PurityComparer(config, Collections.emptyMap());
        MutationBurdenComparer mutationBurdenComparer = new MutationBurdenComparer(config, Collections.emptyMap());

        List<Mismatch> mismatches = Lists.newArrayList();
        assertTrue(purityComparer.processSample(SAMPLE_ID_1, mismatches));
        assertTrue(mutationBurdenComparer.processSample(SAMPLE_ID_1, mismatches));
        assertFalse(purityComparer.processSample(SAMPLE_ID_2, mismatches));

        assertEquals(1, mismatches.size());
        assertEquals(MUTATION_BURDEN, mismatches.get(0).OldItem.category());
        assertEquals(VALUE, mismatches.get(0).Type);
        assertEquals(List.of("TmbPerMb(20.00/2.00)"), mismatches.get(0).DiffValues);
    }
}

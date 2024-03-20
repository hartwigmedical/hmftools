package com.hartwig.hmftools.peach.output;

import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.peach.PeachGenotypeFile;
import com.hartwig.hmftools.peach.HaplotypeAnalysis;
import com.hartwig.hmftools.peach.effect.DrugInfoStore;
import com.hartwig.hmftools.peach.effect.HaplotypeFunctionStore;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.io.IOException;
import java.util.List;
import java.util.Map;

public class BestHaplotypeCombinationsFile
{
    public static void write(@NotNull String filePath, @NotNull Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis,
            @Nullable DrugInfoStore drugInfoStore, @Nullable HaplotypeFunctionStore haplotypeFunctionStore) throws IOException
    {
        List<PeachGenotype> genotypes = PeachGenotypeExtractor.extract(geneToHaplotypeAnalysis, drugInfoStore, haplotypeFunctionStore);
        PeachGenotypeFile.write(filePath, genotypes);
    }
}

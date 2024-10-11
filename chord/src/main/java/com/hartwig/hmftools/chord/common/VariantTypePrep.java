package com.hartwig.hmftools.chord.common;

import java.nio.file.NoSuchFileException;
import java.util.List;

public interface VariantTypePrep<T>
{
    List<T> loadVariants(String sampleId) throws NoSuchFileException;

    List<MutContextCount> countMutationContexts(String sampleId);
}

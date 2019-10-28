package com.hartwig.hmftools.amber;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.amber.ModifiableBaseDepth;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

@Deprecated
public class RefEnricher implements AutoCloseable {

    private final IndexedFastaSequenceFile indexedFastaSequenceFile;

    RefEnricher(@NotNull final String refGenomePath) throws FileNotFoundException {
        indexedFastaSequenceFile = new IndexedFastaSequenceFile(new File(refGenomePath));
    }

    @NotNull
    public ModifiableBaseDepth enrich(@NotNull final ModifiableBaseDepth evidence) {
        try {
            final String baseString =
                    indexedFastaSequenceFile.getSubsequenceAt(evidence.chromosome(), evidence.position(), evidence.position())
                            .getBaseString();
            final BaseDepth.Base base = BaseDepth.Base.valueOf(baseString);
            return evidence.setRef(base);
        } catch (Exception e) {
            return evidence.setRef(BaseDepth.Base.N);
        }
    }

    @Override
    public void close() throws IOException {
        indexedFastaSequenceFile.close();
    }
}

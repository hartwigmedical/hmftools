package com.hartwig.hmftools.amber;

import java.io.File;
import java.io.FileNotFoundException;

import com.hartwig.hmftools.common.amber.ModifiableNormalBAF;
import com.hartwig.hmftools.common.amber.NormalBAF;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class RefEnricher implements AutoCloseable {

    private final IndexedFastaSequenceFile indexedFastaSequenceFile;

    public RefEnricher(@NotNull final String refGenomePath) throws FileNotFoundException {
        indexedFastaSequenceFile = new IndexedFastaSequenceFile(new File(refGenomePath));
    }

    @NotNull
    public ModifiableNormalBAF enrich(@NotNull final ModifiableNormalBAF evidence) {
        try {
            final String baseString =
                    indexedFastaSequenceFile.getSubsequenceAt(evidence.chromosome(), evidence.position(), evidence.position())
                            .getBaseString();
            final NormalBAF.Base base = NormalBAF.Base.valueOf(baseString);
            return evidence.setRef(base);
        } catch (Exception e) {
            return evidence.setRef(NormalBAF.Base.N);
        }
    }

    @Override
    public void close() throws Exception {
        indexedFastaSequenceFile.close();
    }
}

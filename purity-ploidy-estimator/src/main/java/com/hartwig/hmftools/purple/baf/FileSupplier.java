package com.hartwig.hmftools.purple.baf;

import java.io.IOException;
import java.util.function.Supplier;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.purple.baf.TumorBAF;
import com.hartwig.hmftools.common.purple.baf.TumorBAFFile;

class FileSupplier implements Supplier<Multimap<String, TumorBAF>> {

    private final Multimap<String, TumorBAF> bafs;

    FileSupplier(final String filename) throws IOException {
        bafs = TumorBAFFile.read(filename);
    }

    @Override
    public Multimap<String, TumorBAF> get() {
        return bafs;
    }
}

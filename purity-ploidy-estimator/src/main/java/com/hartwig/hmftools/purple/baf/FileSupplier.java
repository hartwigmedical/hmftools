package com.hartwig.hmftools.purple.baf;

import java.util.function.Supplier;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.purple.baf.TumorBAF;

class FileSupplier implements Supplier<Multimap<String, TumorBAF>> {

    private final String filename;

    public FileSupplier(final String filename) {
        this.filename = filename;
    }

    @Override
    public Multimap<String, TumorBAF> get() {
        return null;
    }
}

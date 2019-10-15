package com.hartwig.hmftools.sage.context;

import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.SageVCF;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ContigContext {

    private static final Logger LOGGER = LogManager.getLogger(ContigContext.class);

    private final String contig;
    private final List<Future<List<List<AltContext>>>> futures;

    public ContigContext(final String contig) {
        this.contig = contig;
        futures = Lists.newArrayList();
    }

    public void add(@NotNull final Future<List<List<AltContext>>> futures) {
        this.futures.add(futures);
    }

    public void write(@NotNull final SageVCF vcf) throws ExecutionException, InterruptedException {
        for (Future<List<List<AltContext>>> future : futures) {
            future.get().forEach(vcf::write);
        }

        LOGGER.info("Finalised writing chromosome {} ", contig);
    }

}

package com.hartwig.hmftools.patientdb.dao;

import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.patientdb.Config;

public class SomaticVariantStreamWriter implements Consumer<SomaticVariant> {

    private final String sample;
    private final Timestamp timestamp;
    private final SomaticVariantDAO somaticVariantDAO;
    private final List<SomaticVariant> buffer = new ArrayList<>(Config.DB_BATCH_INSERT_SIZE + 1);

    public SomaticVariantStreamWriter(final SomaticVariantDAO somaticVariantDAO, final String sample) {
        this.somaticVariantDAO = somaticVariantDAO;
        this.sample = sample;
        this.timestamp = new Timestamp(new Date().getTime());
        somaticVariantDAO.deleteSomaticVariantForSample(sample);
    }

    @Override
    public void accept(final SomaticVariant somaticVariant) {
        buffer.add(somaticVariant);
        if (buffer.size() >= Config.DB_BATCH_INSERT_SIZE) {
            writeBuffer();
        }
    }

    public void flush() {
        if (!buffer.isEmpty()) {
            writeBuffer();
        }
    }

    private void writeBuffer() {
        somaticVariantDAO.writeAll(timestamp, sample, buffer);
        buffer.clear();
    }
}

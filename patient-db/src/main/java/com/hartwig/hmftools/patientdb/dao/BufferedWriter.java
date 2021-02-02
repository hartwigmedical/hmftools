package com.hartwig.hmftools.patientdb.dao;

import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.patientdb.Config;

public class BufferedWriter<T> implements Consumer<T>, AutoCloseable {

    private final BufferedWriterConsumer<T> consumer;
    private final Timestamp timestamp;
    private final List<T> buffer = new ArrayList<>(Config.DB_BATCH_INSERT_SIZE + 1);
    private boolean initialised;

    public BufferedWriter(final BufferedWriterConsumer<T> consumer) {
        this.consumer = consumer;
        this.timestamp = new Timestamp(new Date().getTime());
    }

    public void initialise() {
        initialised = true;
        consumer.initialise();
    }

    @Override
    public void accept(final T entry) {
        if (!initialised) {
            initialised = true;
            consumer.initialise();
        }

        buffer.add(entry);
        if (buffer.size() >= Config.DB_BATCH_INSERT_SIZE) {
            writeBuffer();
        }
    }

    private void writeBuffer() {
        consumer.accept(timestamp, buffer);
        buffer.clear();
    }

    @Override
    public void close() {
        if (!buffer.isEmpty()) {
            writeBuffer();
        }
    }
}

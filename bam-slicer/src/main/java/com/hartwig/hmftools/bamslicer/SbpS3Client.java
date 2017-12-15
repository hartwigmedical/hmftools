package com.hartwig.hmftools.bamslicer;

import java.util.concurrent.TimeUnit;

import org.jetbrains.annotations.NotNull;

import okhttp3.ConnectionPool;
import okhttp3.Dispatcher;
import okhttp3.OkHttpClient;

final class SbpS3Client {
    private SbpS3Client() {
    }

    @NotNull
    static OkHttpClient create(final int maxRequests) {
        final Dispatcher requestDispatcher = new Dispatcher();
        requestDispatcher.setMaxRequests(maxRequests);
        requestDispatcher.setMaxRequestsPerHost(maxRequests);
        return new OkHttpClient.Builder().connectionPool(new ConnectionPool(20, 1, TimeUnit.MINUTES))
                .readTimeout(20, TimeUnit.SECONDS)
                .connectTimeout(20, TimeUnit.SECONDS)
                .writeTimeout(20, TimeUnit.SECONDS)
                .dispatcher(requestDispatcher)
                .build();
    }
}

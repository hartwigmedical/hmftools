package com.hartwig.hmftools.bamslicer;

import java.io.IOException;
import java.net.URL;
import java.util.List;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.TimeUnit;

import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.util.concurrent.Futures;
import com.google.common.util.concurrent.ListenableFuture;
import com.google.common.util.concurrent.SettableFuture;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.Chunk;
import htsjdk.samtools.util.BlockCompressedFilePointerUtil;
import okhttp3.Call;
import okhttp3.Callback;
import okhttp3.ConnectionPool;
import okhttp3.Dispatcher;
import okhttp3.Headers;
import okhttp3.OkHttpClient;
import okhttp3.Request;
import okhttp3.Response;
import okhttp3.ResponseBody;

class ChunkHttpBuffer {
    private static final int MAX_REQUESTS = 50;
    @NotNull
    private final LoadingCache<Long, ListenableFuture<byte[]>> chunkBuffer;
    private final int maxSize;
    @NotNull
    private final URL url;
    @NotNull
    private final ConcurrentSkipListMap<Long, Chunk> chunksPerOffset = new ConcurrentSkipListMap<>();
    @NotNull
    private final static OkHttpClient httpClient;

    static {
        final Dispatcher requestDispatcher = new Dispatcher();
        requestDispatcher.setMaxRequests(MAX_REQUESTS);
        requestDispatcher.setMaxRequestsPerHost(MAX_REQUESTS);
        httpClient = new OkHttpClient.Builder().connectionPool(new ConnectionPool(20, 1, TimeUnit.MINUTES))
                .dispatcher(requestDispatcher)
                .build();
    }

    ChunkHttpBuffer(@NotNull final URL url, final int maxSize, @NotNull final List<Chunk> chunks) {
        this.url = url;
        this.maxSize = maxSize;
        for (final Chunk chunk : chunks) {
            final long chunkStart = BlockCompressedFilePointerUtil.getBlockAddress(chunk.getChunkStart());
            chunksPerOffset.put(chunkStart, chunk);
        }
        chunkBuffer = CacheBuilder.newBuilder().maximumSize(maxSize).build(new CacheLoader<Long, ListenableFuture<byte[]>>() {
            @Override
            @NotNull
            public ListenableFuture<byte[]> load(@NotNull final Long offset) throws Exception {
                final Chunk chunkAtOffset = chunksPerOffset.get(offset);
                return getBytesForChunk(chunkAtOffset);
            }

            @Override
            public ListenableFuture<ListenableFuture<byte[]>> reload(@NotNull final Long offset,
                    @NotNull final ListenableFuture<byte[]> oldBytes) {
                return Futures.immediateFuture(oldBytes);
            }
        });
    }

    @NotNull
    Pair<Long, byte[]> getEntryAtPosition(final long position) throws IOException {
        final long chunkOffset = chunksPerOffset.floorEntry(position).getKey();
        final byte[] bytesAtOffset;
        try {
            if (chunkBuffer.getIfPresent(chunkOffset) == null) {
                refillBuffer(chunkOffset);
            }
            bytesAtOffset = chunkBuffer.getUnchecked(chunkOffset).get();
            chunkBuffer.invalidate(chunkOffset);
            if (chunkBuffer.size() < maxSize / 5 + 1) {
                refillBuffer(chunkOffset);
            }
            return Pair.of(chunkOffset, bytesAtOffset);
        } catch (InterruptedException | ExecutionException e) {
            throw new IOException("Could not read entry at position " + position + ". Cause: " + e.getMessage());
        }
    }

    private void refillBuffer(final long chunkOffset) {
        final int fillSize = (int) (maxSize * .75);
        chunksPerOffset.tailMap(chunkOffset, false).keySet().stream().limit(fillSize).forEach(key -> {
            if (chunkBuffer.getIfPresent(key) == null) {
                chunkBuffer.refresh(key);
            }
        });
    }

    @NotNull
    private ListenableFuture<byte[]> getBytesForChunk(@NotNull final Chunk chunk) {
        final long start = BlockCompressedFilePointerUtil.getBlockAddress(chunk.getChunkStart());
        final long end = BlockCompressedFilePointerUtil.getBlockAddress(chunk.getChunkEnd());
        if (start <= end) {
            return readUrlBytes(start, end - start);
        } else {
            return Futures.immediateFailedFuture(new IllegalArgumentException("start offset is greater than end"));
        }
    }

    @NotNull
    private ListenableFuture<byte[]> readUrlBytes(final long offset, final long count) {
        final Headers httpHeaders = new Headers.Builder().add("Range", "bytes=" + offset + "-" + (offset + count - 1)).build();
        final Request request = new Request.Builder().url(url).headers(httpHeaders).build();
        final SettableFuture<byte[]> bytesFuture = SettableFuture.create();
        httpClient.newCall(request).enqueue(retryingCallback(5, bytesFuture));
        return bytesFuture;
    }

    @NotNull
    private Callback retryingCallback(final int retryCount, @NotNull final SettableFuture<byte[]> resultFuture) {
        return new Callback() {
            @Override
            public void onFailure(@NotNull final Call call, @NotNull final IOException e) {
                if (retryCount <= 0) {
                    resultFuture.setException(e);
                } else {
                    call.clone().enqueue(retryingCallback(retryCount - 1, resultFuture));
                }
            }

            @Override
            public void onResponse(@NotNull final Call call, @NotNull final Response response) throws IOException {
                final ResponseBody body = response.body();
                try {
                    if (response.isSuccessful() && body != null) {
                        resultFuture.set(body.bytes());
                    } else {
                        resultFuture.setException(new IOException("Response was not successful or body was null."));
                    }
                } catch (final Exception e) {
                    resultFuture.setException(e);
                } finally {
                    if (body != null) {
                        body.close();
                    }
                }
            }
        };
    }

    @NotNull
    URL url() {
        return url;
    }

    void closeHttpClient() {
        httpClient.dispatcher().executorService().shutdown();
    }
}

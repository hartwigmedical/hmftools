package com.hartwig.hmftools.apiclients.iclusion.http

import okhttp3.ConnectionPool
import okhttp3.Dispatcher
import okhttp3.OkHttpClient
import java.util.concurrent.TimeUnit

private fun requestDispatcher(): Dispatcher {
    val requestDispatcher = Dispatcher()
    requestDispatcher.maxRequests = 100
    requestDispatcher.maxRequestsPerHost = 100
    return requestDispatcher
}

fun httpClient(): OkHttpClient {
    return OkHttpClient.Builder().connectionPool(ConnectionPool(20, 5, TimeUnit.SECONDS))
            .readTimeout(20, TimeUnit.SECONDS)
            .connectTimeout(20, TimeUnit.SECONDS)
            .writeTimeout(20, TimeUnit.SECONDS)
            .dispatcher(requestDispatcher())
            .build()
}

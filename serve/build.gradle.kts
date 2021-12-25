
plugins {
    id("com.hartwig.java-conventions")
}

description = "HMF Tools - SERVE"

dependencies {
    api(project(":hmf-common"))
    implementation(project(":ckb-importer"))
    implementation(project(":iclusion-importer"))
    implementation(project(":vicc-importer"))
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
    testImplementation(libs.junit)
}

// create a jar for tests
testsJar {
}

shadowJar {
    mainClass.set("com.hartwig.hmftools.serve.ServeApplication")
}

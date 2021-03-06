test_that("Create IRIS3 object, expect to be a s4 class IRIS3", {
  data("yan_2013")
  seurat_obj <- Seurat::CreateSeuratObject(yan_2013$expr)
  expect_s4_class(CreateIRIS3Object(seurat_obj), 'IRIS3')
})


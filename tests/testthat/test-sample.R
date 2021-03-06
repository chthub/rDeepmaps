test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("Hello name to be 'Hello: world'", {
  expect_equal(hello_name("world"), 'Hello, world')
})

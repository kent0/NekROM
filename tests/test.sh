# MOR Tests

cd $ROOT_DIR/cases/baf
$ROOT_DIR/bin/linkc

echo 'baf'     > SESSION.NAME
echo `pwd`'/' >> SESSION.NAME

$ROOT_DIR/bin/gnaps baf

case "$TEST" in
    GRAMMIAN_TEST)
        cp $ROOT_DIR/tests/grammain_test.f t.f
        ./nek5000
        ;;

    *)
        echo 'did not specify test...'
        ;;
esac

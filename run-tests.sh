kill $(cat metwork_worker.pid)
rm metwork_worker.pid
bash run-worker.sh test
echo "run unit tests"
./manage.py test --exclude-tag integration --failfast base.tests fragmentation.tests metabolization.tests
echo "restart worker"
kill $(cat metwork_worker.pid)
rm metwork_worker.pid
bash run-worker.sh test
echo "run integration tests"
./manage.py test --tag integration --failfast base.tests
kill $(cat metwork_worker.pid)